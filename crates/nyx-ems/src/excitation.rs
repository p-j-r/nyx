// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Analytical background excitations and the scattered-field formulation of the driven EM study (FR-36, ADR 138).
//!
//! The total field is split as $\mathbf{E} = \mathbf{E}_{\text{inc}} + \mathbf{E}_{\text{scat}}$, where the incident field
//! $\mathbf{E}_{\text{inc}}$ is an exact (or, for the paraxial models, documented approximate) solution of Maxwell's equations in a
//! homogeneous background $(\varepsilon_{\text{bg}}, \mu_{\text{bg}})$. Subtracting the background equation from the total-field
//! curl-curl equation leaves
//! $$\mathbf{A}(\omega)\,\mathbf{E}_{\text{scat}} = \mathbf{b}_{\text{scat}}(\omega),\qquad
//! \langle \mathbf{w}, \mathbf{b}_{\text{scat}}\rangle = \int_{\Omega_{\text{scat}}} \Big[k_0^2\, \mathbf{w}\cdot(\varepsilon_r - \varepsilon_{\text{bg}})\mathbf{E}_{\text{inc}} -
//! (\nabla\times\mathbf{w})\cdot(\mu_r^{-1} - \mu_{\text{bg}}^{-1})(\nabla\times\mathbf{E}_{\text{inc}})\Big]\,dV,$$
//! which vanishes identically wherever the material equals the background (in particular in every PML cell, which is excluded).
//!
//! Time convention $e^{+i\omega t}$ (fields propagate as $e^{-i\mathbf{k}\cdot\mathbf{r}}$, $\mathbf{H} = \tfrac{i}{\omega\mu}\nabla\times\mathbf{E}$),
//! SI units, wafer frame. Excitations:
//! * `PlaneWave`: $\mathbf{E} = A e^{i\phi}\hat{\mathbf{p}}\,e^{-ik\hat{\mathbf{k}}\cdot(\mathbf{r}-\mathbf{r}_0)}$ (exact).
//! * `GaussianBeam` (`complex_source`, default): the regular complex-source/sink beam
//!   $\mathbf{E} \propto \nabla\times\nabla\times(\hat{\mathbf{p}}\, j_0(kR))/k^2$, $R^2 = (\mathbf{r}-\mathbf{r}_f + i b\hat{\mathbf{k}})^2$,
//!   $b = k w_0^2/2$: an exact, singularity-free Maxwell solution whose focal profile is the Gaussian of waist $w_0$.
//!   `paraxial` is the first-order paraxial Gaussian with its longitudinal correction (approximate).
//! * `AstigmaticGaussianBeam`: paraxial elliptical Gaussian with waists $(w_{0x}, w_{0y})$ at axial foci $(z_x, z_y)$ (approximate).
//! * `ThinLensBeam`: vectorial angular-spectrum focusing of an aplanatic lens of numerical aperture NA (a finite Gauss-Legendre x
//!   trapezoid superposition of exact plane waves, with optional Gaussian pupil filling).
//! * `Dipole`: Hertzian electric or magnetic point dipole radiating into the background (exact; must lie outside the contrast region).

use crate::cavity::{basis3_keyed, safe_curved_jacobian, Ned3, CURVED_QUAD_DEGREE};
pub use crate::ibc::EPS0;
use crate::ibc::MU0;
use crate::modes::C0;
use femcore::fem::{geometry, quad_rule, quad_rule_collapsed, CurvedTet, Geometry};
use femcore::linalg::csr::{c, C};
use nyx_physics::Discretisation;
use std::f64::consts::PI;

type C3 = [[C; 3]; 3];

/// Field model of a `GaussianBeam`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BeamModel {
    /// Exact complex-source/sink beam (default).
    ComplexSource,
    /// First-order paraxial Gaussian with longitudinal correction.
    Paraxial,
}

/// One analytical incident field.
#[derive(Clone, Debug, PartialEq)]
pub enum Excitation {
    PlaneWave {
        k_hat: [f64; 3],
        pol: [C; 3],
        amplitude: f64,
        phase: f64,
        origin: [f64; 3],
    },
    GaussianBeam {
        focus: [f64; 3],
        w0: f64,
        k_hat: [f64; 3],
        pol: [C; 3],
        amplitude: f64,
        phase: f64,
        model: BeamModel,
    },
    AstigmaticGaussianBeam {
        focus: [f64; 3],
        w0_x: f64,
        w0_y: f64,
        focus_x: f64,
        focus_y: f64,
        k_hat: [f64; 3],
        x_axis: [f64; 3],
        pol: [C; 3],
        amplitude: f64,
        phase: f64,
    },
    ThinLensBeam {
        focus: [f64; 3],
        na: f64,
        k_hat: [f64; 3],
        x_axis: [f64; 3],
        pol: [C; 3],
        amplitude: f64,
        phase: f64,
        filling: f64,
        n_theta: usize,
        n_phi: usize,
    },
    Dipole {
        position: [f64; 3],
        moment: [C; 3],
        magnetic: bool,
    },
}

/// One scattered-field excitation of a driven study: the incident field, the homogeneous background it solves, and the contrast
/// region $\Omega_{\text{scat}}$ (domain labels; empty = every non-PML cell whose material differs from the background).
#[derive(Clone, Debug, PartialEq)]
pub struct ExcitationSpec {
    pub name: String,
    pub excitation: Excitation,
    pub domains: Vec<i32>,
    pub eps_bg: C,
    pub mu_bg: f64,
    /// Also emit the scattered fields `Esx..Esz`, `Hsx..Hsz` next to the total fields.
    pub keep_scattered: bool,
}

#[inline]
fn dot_r(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn cross_r(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
fn cross_c(a: [C; 3], b: [C; 3]) -> [C; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
#[inline]
fn dot_c(a: [C; 3], b: [C; 3]) -> C {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn norm_r(a: [f64; 3]) -> f64 {
    dot_r(a, a).sqrt()
}
#[inline]
fn unit(a: [f64; 3], what: &str) -> Result<[f64; 3], String> {
    let n = norm_r(a);
    if !(n > 0.0 && n.is_finite()) {
        return Err(format!("{what} must be a non-zero finite vector"));
    }
    Ok([a[0] / n, a[1] / n, a[2] / n])
}
#[inline]
fn rc(a: [f64; 3]) -> [C; 3] {
    [c(a[0], 0.0), c(a[1], 0.0), c(a[2], 0.0)]
}
#[inline]
fn scale_c(s: C, a: [C; 3]) -> [C; 3] {
    [s * a[0], s * a[1], s * a[2]]
}
#[inline]
fn add_c(a: [C; 3], b: [C; 3]) -> [C; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Transverse polarization: removes the component along `k_hat` and normalizes to unit Hermitian norm.
fn transverse_pol(pol: [C; 3], k_hat: [f64; 3]) -> Result<[C; 3], String> {
    let kc = rc(k_hat);
    let along = dot_c(pol, kc);
    let p = [pol[0] - along * kc[0], pol[1] - along * kc[1], pol[2] - along * kc[2]];
    let n = (p[0].norm_sqr() + p[1].norm_sqr() + p[2].norm_sqr()).sqrt();
    if !(n > 1e-12 && n.is_finite()) {
        return Err("excitation polarization has no component transverse to the propagation direction".into());
    }
    Ok(scale_c(c(1.0 / n, 0.0), p))
}

/// Orthonormal transverse frame `(u, v)` with `u` along the projection of `x_axis` and `v = k x u`.
fn transverse_frame(k_hat: [f64; 3], x_axis: [f64; 3]) -> Result<([f64; 3], [f64; 3]), String> {
    let mut u = [
        x_axis[0] - dot_r(x_axis, k_hat) * k_hat[0],
        x_axis[1] - dot_r(x_axis, k_hat) * k_hat[1],
        x_axis[2] - dot_r(x_axis, k_hat) * k_hat[2],
    ];
    if norm_r(u) < 1e-9 {
        // x_axis parallel to k: fall back to the global axis least aligned with k
        let alt = if k_hat[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        u = [
            alt[0] - dot_r(alt, k_hat) * k_hat[0],
            alt[1] - dot_r(alt, k_hat) * k_hat[1],
            alt[2] - dot_r(alt, k_hat) * k_hat[2],
        ];
    }
    let u = unit(u, "transverse axis")?;
    let v = cross_r(k_hat, u);
    Ok((u, v))
}

/// Gauss-Legendre nodes and weights on `[-1, 1]` (Newton iteration on the Legendre recurrence).
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let n = n.max(1);
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, 0.0);
            for j in 0..n {
                let p2 = p1;
                p1 = p0;
                p0 = ((2.0 * j as f64 + 1.0) * z * p1 - j as f64 * p2) / (j as f64 + 1.0);
            }
            dp = n as f64 * (z * p0 - p1) / (z * z - 1.0);
            let dz = p0 / dp;
            z -= dz;
            if dz.abs() < 1e-15 {
                break;
            }
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
        w[n - 1 - i] = w[i];
    }
    (x, w)
}

impl Excitation {
    /// Checks the parameters (finite, positive sizes, NA below the background index, ...).
    pub fn validate(&self, n_bg: f64) -> Result<(), String> {
        let finite3 = |v: &[f64; 3]| v.iter().all(|x| x.is_finite());
        let finite_c3 = |v: &[C; 3]| v.iter().all(|z| z.re.is_finite() && z.im.is_finite());
        match self {
            Excitation::PlaneWave { k_hat, pol, amplitude, phase, origin } => {
                if !(finite3(k_hat) && finite_c3(pol) && amplitude.is_finite() && phase.is_finite() && finite3(origin))
                {
                    return Err("plane_wave parameters must be finite".into());
                }
                transverse_pol(*pol, unit(*k_hat, "plane_wave direction")?).map(|_| ())
            }
            Excitation::GaussianBeam { focus, w0, k_hat, pol, amplitude, phase, .. } => {
                if !(finite3(focus) && finite3(k_hat) && finite_c3(pol) && amplitude.is_finite() && phase.is_finite()) {
                    return Err("gaussian_beam parameters must be finite".into());
                }
                if !(*w0 > 0.0 && w0.is_finite()) {
                    return Err(format!("gaussian_beam waist w0 must be positive (got {w0})"));
                }
                transverse_pol(*pol, unit(*k_hat, "gaussian_beam direction")?).map(|_| ())
            }
            Excitation::AstigmaticGaussianBeam {
                focus,
                w0_x,
                w0_y,
                focus_x,
                focus_y,
                k_hat,
                x_axis,
                pol,
                amplitude,
                phase,
            } => {
                if !(finite3(focus)
                    && finite3(k_hat)
                    && finite3(x_axis)
                    && finite_c3(pol)
                    && amplitude.is_finite()
                    && phase.is_finite()
                    && focus_x.is_finite()
                    && focus_y.is_finite())
                {
                    return Err("astigmatic_gaussian_beam parameters must be finite".into());
                }
                if !(*w0_x > 0.0 && *w0_y > 0.0 && w0_x.is_finite() && w0_y.is_finite()) {
                    return Err("astigmatic_gaussian_beam waists w0_x, w0_y must be positive".into());
                }
                transverse_pol(*pol, unit(*k_hat, "astigmatic_gaussian_beam direction")?).map(|_| ())
            }
            Excitation::ThinLensBeam { focus, na, k_hat, x_axis, pol, amplitude, phase, filling, n_theta, n_phi } => {
                if !(finite3(focus)
                    && finite3(k_hat)
                    && finite3(x_axis)
                    && finite_c3(pol)
                    && amplitude.is_finite()
                    && phase.is_finite()
                    && filling.is_finite())
                {
                    return Err("thin_lens_beam parameters must be finite".into());
                }
                if !(*na > 0.0 && *na < n_bg) {
                    return Err(format!("thin_lens_beam requires 0 < NA < n_background = {n_bg:.6} (got {na})"));
                }
                if *n_theta < 2 || *n_phi < 4 {
                    return Err("thin_lens_beam requires n_theta >= 2 and n_phi >= 4".into());
                }
                transverse_pol(*pol, unit(*k_hat, "thin_lens_beam direction")?).map(|_| ())
            }
            Excitation::Dipole { position, moment, .. } => {
                if !(finite3(position) && finite_c3(moment)) {
                    return Err("dipole parameters must be finite".into());
                }
                if moment.iter().all(|z| z.norm() == 0.0) {
                    return Err("dipole moment must be non-zero".into());
                }
                Ok(())
            }
        }
    }

    /// Whether the incident field is a plane wave (cross-sections are then normalized by its intensity).
    pub fn is_plane_wave(&self) -> bool {
        matches!(self, Excitation::PlaneWave { .. })
    }

    /// Prepares the field at angular frequency `omega` in the background `(eps_bg, mu_bg)`.
    pub fn prepare(&self, omega: f64, eps_bg: C, mu_bg: f64) -> Result<PreparedExcitation, String> {
        let n_c = (eps_bg * mu_bg).sqrt();
        // e^{+i w t}: a lossy background has Im(n) <= 0
        let n_c = if n_c.im > 0.0 { -n_c } else { n_c };
        let k = n_c * (omega / C0);
        self.validate(n_c.re)?;
        let kind = match self {
            Excitation::PlaneWave { k_hat, pol, amplitude, phase, origin } => {
                let kh = unit(*k_hat, "plane_wave direction")?;
                let p = transverse_pol(*pol, kh)?;
                let a = c(amplitude * phase.cos(), amplitude * phase.sin());
                Prepared::Plane { k_hat: kh, e0: scale_c(a, p), origin: *origin }
            }
            Excitation::GaussianBeam { focus, w0, k_hat, pol, amplitude, phase, model } => {
                let kh = unit(*k_hat, "gaussian_beam direction")?;
                let p = transverse_pol(*pol, kh)?;
                match model {
                    BeamModel::ComplexSource => {
                        let b = 0.5 * k.re * w0 * w0;
                        Prepared::Csp { focus: *focus, k_hat: kh, b, kappa: k.re * b, p, norm: c(1.0, 0.0) }
                    }
                    BeamModel::Paraxial => {
                        let (u, v) = transverse_frame(kh, [1.0, 0.0, 0.0])?;
                        Prepared::Paraxial {
                            focus: *focus,
                            k_hat: kh,
                            u,
                            v,
                            zr: [0.5 * k.re * w0 * w0, 0.5 * k.re * w0 * w0],
                            zf: [0.0, 0.0],
                            pu: dot_c(p, rc(u)),
                            pv: dot_c(p, rc(v)),
                            norm: c(1.0, 0.0),
                        }
                    }
                }
                .with_target(*amplitude, *phase, p, k, omega, *focus)?
            }
            Excitation::AstigmaticGaussianBeam {
                focus,
                w0_x,
                w0_y,
                focus_x,
                focus_y,
                k_hat,
                x_axis,
                pol,
                amplitude,
                phase,
            } => {
                let kh = unit(*k_hat, "astigmatic_gaussian_beam direction")?;
                let p = transverse_pol(*pol, kh)?;
                let (u, v) = transverse_frame(kh, *x_axis)?;
                Prepared::Paraxial {
                    focus: *focus,
                    k_hat: kh,
                    u,
                    v,
                    zr: [0.5 * k.re * w0_x * w0_x, 0.5 * k.re * w0_y * w0_y],
                    zf: [*focus_x, *focus_y],
                    pu: dot_c(p, rc(u)),
                    pv: dot_c(p, rc(v)),
                    norm: c(1.0, 0.0),
                }
                .with_target(*amplitude, *phase, p, k, omega, *focus)?
            }
            Excitation::ThinLensBeam { focus, na, k_hat, x_axis, pol, amplitude, phase, filling, n_theta, n_phi } => {
                let kh = unit(*k_hat, "thin_lens_beam direction")?;
                let p = transverse_pol(*pol, kh)?;
                let (u, v) = transverse_frame(kh, *x_axis)?;
                let (pu, pv) = (dot_c(p, rc(u)), dot_c(p, rc(v)));
                let theta_max = (na / n_c.re).clamp(0.0, 1.0).asin();
                let (xg, wg) = gauss_legendre(*n_theta);
                let mut waves = Vec::with_capacity(n_theta * n_phi);
                for (xi, wi) in xg.iter().zip(&wg) {
                    let th = 0.5 * theta_max * (xi + 1.0);
                    let w_th = 0.5 * theta_max * wi;
                    let (st, ct) = th.sin_cos();
                    let apod_fill = if *filling > 0.0 {
                        (-(st / theta_max.sin()).powi(2) / (filling * filling)).exp()
                    } else {
                        1.0
                    };
                    let weight = w_th * st * ct.sqrt() * apod_fill * (2.0 * PI / *n_phi as f64);
                    for j in 0..*n_phi {
                        let ph = 2.0 * PI * j as f64 / *n_phi as f64;
                        let (sp, cp) = ph.sin_cos();
                        let es_t = [-sp, cp]; // in (u, v)
                        let ep_in = [cp, sp];
                        let a_s = pu * es_t[0] + pv * es_t[1];
                        let a_p = pu * ep_in[0] + pv * ep_in[1];
                        let mut e = [c(0.0, 0.0); 3];
                        for k3 in 0..3 {
                            let es = es_t[0] * u[k3] + es_t[1] * v[k3];
                            let ep = ct * cp * u[k3] + ct * sp * v[k3] - st * kh[k3];
                            e[k3] = (a_s * es + a_p * ep) * weight;
                        }
                        let s_hat = [
                            st * cp * u[0] + st * sp * v[0] + ct * kh[0],
                            st * cp * u[1] + st * sp * v[1] + ct * kh[1],
                            st * cp * u[2] + st * sp * v[2] + ct * kh[2],
                        ];
                        waves.push((e, s_hat));
                    }
                }
                Prepared::Lens { focus: *focus, waves, norm: c(1.0, 0.0) }
                    .with_target(*amplitude, *phase, p, k, omega, *focus)?
            }
            Excitation::Dipole { position, moment, magnetic } => {
                Prepared::Dipole { position: *position, moment: *moment, magnetic: *magnetic }
            }
        };
        Ok(PreparedExcitation { kind, k, omega, eps_abs: eps_bg * EPS0, mu_abs: mu_bg * MU0 })
    }
}

/// Internal prepared representation.
#[derive(Clone, Debug)]
enum Prepared {
    Plane {
        k_hat: [f64; 3],
        e0: [C; 3],
        origin: [f64; 3],
    },
    Csp {
        focus: [f64; 3],
        k_hat: [f64; 3],
        b: f64,
        kappa: f64,
        p: [C; 3],
        norm: C,
    },
    Paraxial {
        focus: [f64; 3],
        k_hat: [f64; 3],
        u: [f64; 3],
        v: [f64; 3],
        zr: [f64; 2],
        zf: [f64; 2],
        pu: C,
        pv: C,
        norm: C,
    },
    Lens {
        focus: [f64; 3],
        waves: Vec<([C; 3], [f64; 3])>,
        norm: C,
    },
    Dipole {
        position: [f64; 3],
        moment: [C; 3],
        magnetic: bool,
    },
}

impl Prepared {
    /// Scales the raw field so that `conj(p) . E(focus) = amplitude e^{i phase}`.
    fn with_target(
        mut self,
        amplitude: f64,
        phase: f64,
        p: [C; 3],
        k: C,
        omega: f64,
        focus: [f64; 3],
    ) -> Result<Prepared, String> {
        let probe = PreparedExcitation { kind: self.clone(), k, omega, eps_abs: c(EPS0, 0.0), mu_abs: MU0 };
        let (e, _) = probe.field(focus);
        let proj = p[0].conj() * e[0] + p[1].conj() * e[1] + p[2].conj() * e[2];
        if !(proj.norm() > 1e-300 && proj.norm().is_finite()) {
            return Err("beam excitation has a vanishing field at its focus: cannot normalize its amplitude".into());
        }
        let target = c(amplitude * phase.cos(), amplitude * phase.sin()) / proj;
        match &mut self {
            Prepared::Csp { norm, .. } | Prepared::Paraxial { norm, .. } | Prepared::Lens { norm, .. } => {
                *norm = target
            }
            _ => {}
        }
        Ok(self)
    }
}

/// An excitation evaluated at one frequency: `field(x) -> (E, curl E)` in V/m and V/m^2.
#[derive(Clone, Debug)]
pub struct PreparedExcitation {
    kind: Prepared,
    /// Background wavenumber (1/m, `Im <= 0` when lossy).
    pub k: C,
    pub omega: f64,
    /// Absolute background permittivity `eps0 eps_bg` (F/m).
    pub eps_abs: C,
    /// Absolute background permeability `mu0 mu_bg` (H/m).
    pub mu_abs: f64,
}

/// Scaled spherical Bessel ratios `(j0(z), j1(z)/z, j2(z)/z^2) * exp(-kappa)` (entire functions of `z^2`).
fn sph_bessel_ratios_scaled(z: C, kappa: f64) -> (C, C, C) {
    if z.norm() < 1.0 {
        let s = (-kappa).exp();
        let z2h = -(z * z) * 0.5;
        let series = |n: usize| -> C {
            // j_n(z)/z^n = sum_m (-z^2/2)^m / (m! (2m+2n+1)!!)
            let mut df = 1.0;
            for j in (1..=(2 * n + 1)).step_by(2) {
                df *= j as f64;
            }
            let mut term = c(1.0 / df, 0.0);
            let mut sum = term;
            for m in 1..30 {
                term = term * z2h / ((m as f64) * (2 * m + 2 * n + 1) as f64);
                sum += term;
                if term.norm() < 1e-18 * sum.norm() {
                    break;
                }
            }
            sum * s
        };
        return (series(0), series(1), series(2));
    }
    let iz = c(0.0, 1.0) * z;
    let ep = (iz - kappa).exp();
    let em = (-iz - kappa).exp();
    let sin_s = (ep - em) / c(0.0, 2.0);
    let cos_s = (ep + em) * 0.5;
    let z2 = z * z;
    let j0 = sin_s / z;
    let j1z = (sin_s - z * cos_s) / (z2 * z);
    let j2z2 = ((c(3.0, 0.0) - z2) * sin_s - z * cos_s * 3.0) / (z2 * z2 * z);
    (j0, j1z, j2z2)
}

impl PreparedExcitation {
    /// Incident electric field and its curl at the physical point `x`.
    pub fn field(&self, x: [f64; 3]) -> ([C; 3], [C; 3]) {
        let ci = c(0.0, 1.0);
        let k = self.k;
        match &self.kind {
            Prepared::Plane { k_hat, e0, origin } => {
                let s = (x[0] - origin[0]) * k_hat[0] + (x[1] - origin[1]) * k_hat[1] + (x[2] - origin[2]) * k_hat[2];
                let ph = (-ci * k * s).exp();
                let e = scale_c(ph, *e0);
                let curl = scale_c(-ci * k, cross_c(rc(*k_hat), e));
                (e, curl)
            }
            Prepared::Csp { focus, k_hat, b, kappa, p, norm } => {
                let d = [
                    c(x[0] - focus[0], b * k_hat[0]),
                    c(x[1] - focus[1], b * k_hat[1]),
                    c(x[2] - focus[2], b * k_hat[2]),
                ];
                let r2 = dot_c(d, d);
                let r = r2.sqrt();
                let z = k * r;
                let (j0, j1z, j2z2) = sph_bessel_ratios_scaled(z, *kappa);
                let dp = dot_c(d, *p);
                let k2 = k * k;
                let e = add_c(scale_c(*norm * (j0 - j1z), *p), scale_c(*norm * k2 * j2z2 * dp, d));
                let curl = scale_c(-*norm * k2 * j1z, cross_c(d, *p));
                (e, curl)
            }
            Prepared::Paraxial { .. } => {
                let e = self.paraxial_e(x);
                // central-difference curl (step lambda/2000): O(h^2 k^2) ~ 1e-5 relative
                let h = (2.0 * PI / k.norm().max(1e-300)) / 2000.0;
                let mut grad = [[c(0.0, 0.0); 3]; 3]; // grad[a][comp] = d E_comp / d x_a
                for a in 0..3 {
                    let mut xp = x;
                    let mut xm = x;
                    xp[a] += h;
                    xm[a] -= h;
                    let (ep, em) = (self.paraxial_e(xp), self.paraxial_e(xm));
                    for comp in 0..3 {
                        grad[a][comp] = (ep[comp] - em[comp]) / (2.0 * h);
                    }
                }
                let curl = [grad[1][2] - grad[2][1], grad[2][0] - grad[0][2], grad[0][1] - grad[1][0]];
                (e, curl)
            }
            Prepared::Lens { focus, waves, norm } => {
                let mut e = [c(0.0, 0.0); 3];
                let mut curl = [c(0.0, 0.0); 3];
                let r = [x[0] - focus[0], x[1] - focus[1], x[2] - focus[2]];
                for (amp, s_hat) in waves {
                    let ph = (-ci * k * dot_r(*s_hat, r)).exp() * *norm;
                    let ew = scale_c(ph, *amp);
                    e = add_c(e, ew);
                    curl = add_c(curl, scale_c(-ci * k, cross_c(rc(*s_hat), ew)));
                }
                (e, curl)
            }
            Prepared::Dipole { position, moment, magnetic } => {
                let d = [x[0] - position[0], x[1] - position[1], x[2] - position[2]];
                let r = norm_r(d).max(1e-300);
                let n = rc([d[0] / r, d[1] / r, d[2] / r]);
                let rr = c(r, 0.0);
                let f = (-ci * k * rr).exp() / (4.0 * PI * r);
                let k2 = k * k;
                let np = dot_c(n, *moment);
                // E = curl curl (m f) pattern: (k^2 - i k/R - 1/R^2) m + (-k^2 + 3 i k/R + 3/R^2)(n.m) n, times f
                let a1 = k2 - ci * k / rr - c(1.0 / (r * r), 0.0);
                let a2 = -k2 + ci * k * 3.0 / rr + c(3.0 / (r * r), 0.0);
                let cc = add_c(scale_c(a1 * f, *moment), scale_c(a2 * f * np, n));
                // grad f x m = f'(R) n x m, f' = -(i k + 1/R) f
                let fp = -(ci * k + c(1.0 / r, 0.0)) * f;
                let gxm = scale_c(fp, cross_c(n, *moment));
                if *magnetic {
                    let pref = -ci * self.omega * self.mu_abs;
                    (scale_c(pref, gxm), scale_c(pref, cc))
                } else {
                    let inv_eps = c(1.0, 0.0) / self.eps_abs;
                    (scale_c(inv_eps, cc), scale_c(inv_eps * k2, gxm))
                }
            }
        }
    }

    fn paraxial_e(&self, x: [f64; 3]) -> [C; 3] {
        let Prepared::Paraxial { focus, k_hat, u, v, zr, zf, pu, pv, norm } = &self.kind else {
            return [c(0.0, 0.0); 3];
        };
        let ci = c(0.0, 1.0);
        let k = self.k;
        let r = [x[0] - focus[0], x[1] - focus[1], x[2] - focus[2]];
        let (xu, xv, zeta) = (dot_r(r, *u), dot_r(r, *v), dot_r(r, *k_hat));
        let qx = c(zeta - zf[0], zr[0]);
        let qy = c(zeta - zf[1], zr[1]);
        let q0x = c(0.0, zr[0]);
        let q0y = c(0.0, zr[1]);
        let env = (q0x / qx).sqrt() * (q0y / qy).sqrt();
        let ph = (-ci * k * zeta - ci * k * xu * xu / (qx * 2.0) - ci * k * xv * xv / (qy * 2.0)).exp();
        let s = *norm * env * ph;
        let eu = *pu * s;
        let ev = *pv * s;
        let ez = -(eu * xu / qx) - (ev * xv / qy);
        let mut e = [c(0.0, 0.0); 3];
        for a in 0..3 {
            e[a] = eu * u[a] + ev * v[a] + ez * k_hat[a];
        }
        e
    }

    /// Incident magnetic field `H = i/(omega mu_bg) curl E`.
    pub fn h_field(&self, curl_e: [C; 3]) -> [C; 3] {
        scale_c(c(0.0, 1.0 / (self.omega * self.mu_abs)), curl_e)
    }

    /// Plane-wave intensity `|E0|^2 / (2 eta_bg)` (W/m^2), `None` for other excitations.
    pub fn plane_wave_intensity(&self) -> Option<f64> {
        let Prepared::Plane { e0, .. } = &self.kind else { return None };
        let e2 = e0[0].norm_sqr() + e0[1].norm_sqr() + e0[2].norm_sqr();
        // eta = sqrt(mu/eps); intensity = |E|^2 Re(1/eta*) / 2
        let eta = (c(self.mu_abs, 0.0) / self.eps_abs).sqrt();
        Some(0.5 * e2 * (c(1.0, 0.0) / eta.conj()).re)
    }

    /// Power radiated by a point dipole into the lossless background (W), `None` for other excitations.
    pub fn dipole_power(&self) -> Option<f64> {
        let Prepared::Dipole { moment, magnetic, .. } = &self.kind else { return None };
        let m2 = moment[0].norm_sqr() + moment[1].norm_sqr() + moment[2].norm_sqr();
        let k3 = self.k.re.powi(3);
        Some(if *magnetic {
            self.omega * self.mu_abs * k3 * m2 / (12.0 * PI)
        } else {
            self.omega * k3 * m2 / (12.0 * PI * self.eps_abs.re)
        })
    }

    /// Dipole position (to reject contrast cells that contain it).
    pub fn dipole_position(&self) -> Option<[f64; 3]> {
        match &self.kind {
            Prepared::Dipole { position, .. } => Some(*position),
            _ => None,
        }
    }
}

// ------------------------------------------------------------------ cell quadrature over the contrast region

/// One quadrature point of a tetrahedron: barycentric coordinates, weight (with the Jacobian), physical position and gradients.
pub(crate) struct CellQp {
    pub lam: [f64; 4],
    pub wq: f64,
    pub xq: [f64; 3],
    pub g: [[f64; 3]; 4],
}

/// Quadrature points of cell `ci` (degree-4 Keast rule on straight cells, the curved rule of ADR 82 on curved cells).
pub(crate) fn cell_quadrature(mesh: &Discretisation, geo: &Geometry, ci: usize) -> Vec<CellQp> {
    let mut gs = [[0.0; 3]; 4];
    for b in 0..4 {
        for k in 0..3 {
            gs[b][k] = geo.g[ci * 12 + b * 3 + k];
        }
    }
    let mut out = Vec::new();
    if mesh.cell_is_curved(ci) {
        let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
        let cq = quad_rule_collapsed(CURVED_QUAD_DEGREE);
        for q in 0..cq.w.len() {
            let lam = [cq.lam[q][0], cq.lam[q][1], cq.lam[q][2], cq.lam[q][3]];
            let (_, det) = ct.jacobian(&lam);
            let g = safe_curved_jacobian(&ct, &lam, &gs);
            out.push(CellQp { lam, wq: cq.w[q] * det.abs() / 6.0, xq: ct.position(&lam), g });
        }
    } else {
        let vs = mesh.cell(ci);
        let pts: Vec<&[f64]> = vs.iter().map(|&v| mesh.point(v as usize)).collect();
        let quad = quad_rule(3, 6);
        for q in 0..quad.w.len() {
            let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
            let mut xq = [0.0; 3];
            for b in 0..4 {
                for k in 0..3 {
                    xq[k] += lam[b] * pts[b][k];
                }
            }
            out.push(CellQp { lam, wq: quad.w[q] * geo.vol[ci], xq, g: gs });
        }
    }
    out
}

/// Material contrast of one cell against the background: `(eps - eps_bg I, mu_inv - mu_bg^-1 I)`.
#[inline]
pub(crate) fn contrast(eps: &C3, mu_inv: &C3, eps_bg: C, mu_bg: f64) -> (C3, C3, bool) {
    let mut de = *eps;
    let mut dm = *mu_inv;
    for a in 0..3 {
        de[a][a] -= eps_bg;
        dm[a][a] -= c(1.0 / mu_bg, 0.0);
    }
    let scale = eps_bg.norm().max(1.0);
    let nz = de.iter().flatten().chain(dm.iter().flatten()).any(|z| z.norm() > 1e-12 * scale);
    (de, dm, nz)
}

#[inline]
fn c3_apply(m: &C3, v: [C; 3]) -> [C; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Cells of the contrast region of `spec`: listed domains (or all), non-PML, with a non-zero material contrast.
pub(crate) fn contrast_cells(
    mesh: &Discretisation,
    spec: &ExcitationSpec,
    eps: &[C3],
    mu_inv: &[C3],
    pml_cells: &[bool],
) -> Vec<usize> {
    (0..mesh.n_cells())
        .filter(|&ci| {
            !pml_cells.get(ci).copied().unwrap_or(false)
                && (spec.domains.is_empty() || spec.domains.contains(&mesh.cell_labels[ci]))
                && contrast(&eps[ci], &mu_inv[ci], spec.eps_bg, spec.mu_bg).2
        })
        .collect()
}

/// Whether the physical point `p` lies inside the straight tetrahedron `ci` (curved cells use their vertex hull).
fn point_in_cell(mesh: &Discretisation, ci: usize, p: [f64; 3]) -> bool {
    let vs = mesh.cell(ci);
    let a: Vec<[f64; 3]> = vs
        .iter()
        .map(|&v| {
            let q = mesh.point(v as usize);
            [q[0], q[1], q[2]]
        })
        .collect();
    let sub = |x: [f64; 3], y: [f64; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let vol =
        |p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]| dot_r(sub(p1, p0), cross_r(sub(p2, p0), sub(p3, p0)));
    let v0 = vol(a[0], a[1], a[2], a[3]);
    if v0.abs() < 1e-300 {
        return false;
    }
    let l = [
        vol(p, a[1], a[2], a[3]) / v0,
        vol(a[0], p, a[2], a[3]) / v0,
        vol(a[0], a[1], p, a[3]) / v0,
        vol(a[0], a[1], a[2], p) / v0,
    ];
    l.iter().all(|&x| x >= -1e-9)
}

/// Fixed cell-chunk size of the deterministic scattered-field loops (same order for every `NYX_THREADS`).
const EXC_CHUNK: usize = 1024;

/// Assembles the full-DOF scattered-field load vector $\mathbf{b}_{\text{scat}}(\omega)$ (length `ned.n_dofs()`).
pub fn assemble_scattered_rhs(
    mesh: &Discretisation,
    ned: &Ned3,
    spec: &ExcitationSpec,
    prep: &PreparedExcitation,
    eps: &[C3],
    mu_inv: &[C3],
    pml_cells: &[bool],
) -> Result<Vec<C>, String> {
    let geo = geometry(mesh);
    let cells = contrast_cells(mesh, spec, eps, mu_inv, pml_cells);
    if let Some(pos) = prep.dipole_position() {
        if let Some(&ci) = cells.iter().find(|&&ci| point_in_cell(mesh, ci, pos)) {
            return Err(format!(
                "excitation '{}': the dipole lies inside contrast cell {ci} (domain {}); place it in the background or restrict 'on' to the scatterer domains",
                spec.name, mesh.cell_labels[ci]
            ));
        }
    }
    let k0 = prep.omega / C0;
    let k02 = k0 * k0;
    let chunks: Vec<&[usize]> = cells.chunks(EXC_CHUNK).collect();
    let parts: Vec<Vec<(usize, C)>> =
        femcore::par::chunk_map(chunks.len(), chunks.len() > 1 && femcore::par::threads() > 1, |kc| {
            let mut out = Vec::new();
            let (mut w, mut cu) = (Vec::new(), Vec::new());
            for &ci in chunks[kc] {
                let (de, dm, _) = contrast(&eps[ci], &mu_inv[ci], spec.eps_bg, spec.mu_bg);
                let gid = ned.key_pairs(mesh, ci);
                let nb = ned.nb();
                let mut fe = vec![c(0.0, 0.0); nb];
                for qp in cell_quadrature(mesh, &geo, ci) {
                    basis3_keyed(ned.order, &qp.lam, &qp.g, &gid, &mut w, &mut cu);
                    let (e_inc, curl_inc) = prep.field(qp.xq);
                    let s_e = c3_apply(&de, e_inc);
                    let s_c = c3_apply(&dm, curl_inc);
                    for s in 0..nb {
                        let we = w[s][0] * s_e[0] + w[s][1] * s_e[1] + w[s][2] * s_e[2];
                        let wc = cu[s][0] * s_c[0] + cu[s][1] * s_c[1] + cu[s][2] * s_c[2];
                        fe[s] += (we * k02 - wc) * qp.wq;
                    }
                }
                for s in 0..nb {
                    out.push((ned.cell_dofs[ci][s], fe[s]));
                }
            }
            out
        });
    let mut b = vec![c(0.0, 0.0); ned.n_dofs()];
    for part in parts {
        for (d, v) in part {
            b[d] += v;
        }
    }
    if b.iter().any(|z| !(z.re.is_finite() && z.im.is_finite())) {
        return Err(format!("excitation '{}': the scattered-field load vector is not finite", spec.name));
    }
    Ok(b)
}

/// Volume powers of a solved excitation: `(P_abs, P_ext)` in W.
///
/// * $P_{\text{abs}} = -\tfrac{\omega\varepsilon_0}{2}\int \operatorname{Im}(\mathbf{E}^*\cdot\varepsilon_r\mathbf{E})\,dV$ over every
///   non-PML cell (total field), plus the magnetic analogue for complex `mu_inv`.
/// * $P_{\text{ext}} = -\tfrac{\omega\varepsilon_0}{2}\operatorname{Im}\int_{\Omega_{\text{scat}}} \mathbf{E}_{\text{inc}}^*\cdot\Delta\varepsilon\,\mathbf{E}\,dV -
///   \tfrac{\omega\mu_0}{2}\operatorname{Im}\int \mathbf{H}_{\text{inc}}^*\cdot\Delta\mu\,\mathbf{H}\,dV$ (work of the incident field on
///   the induced polarization currents; the optical theorem gives $P_{\text{sca}} = P_{\text{ext}} - P_{\text{abs}}$).
pub fn excitation_volume_powers(
    mesh: &Discretisation,
    ned: &Ned3,
    spec: &ExcitationSpec,
    prep: &PreparedExcitation,
    eps: &[C3],
    mu_inv: &[C3],
    pml_cells: &[bool],
    x_scat_full: &[C],
) -> (f64, f64) {
    let geo = geometry(mesh);
    let omega = prep.omega;
    let contrast_set = contrast_cells(mesh, spec, eps, mu_inv, pml_cells);
    let mut is_contrast = vec![false; mesh.n_cells()];
    for &ci in &contrast_set {
        is_contrast[ci] = true;
    }
    let lossy = |ci: usize| -> bool {
        eps[ci].iter().flatten().any(|z| z.im.abs() > 0.0) || mu_inv[ci].iter().flatten().any(|z| z.im.abs() > 0.0)
    };
    let cells: Vec<usize> = (0..mesh.n_cells())
        .filter(|&ci| !pml_cells.get(ci).copied().unwrap_or(false) && (is_contrast[ci] || lossy(ci)))
        .collect();
    let chunks: Vec<&[usize]> = cells.chunks(EXC_CHUNK).collect();
    let parts: Vec<(f64, f64)> =
        femcore::par::chunk_map(chunks.len(), chunks.len() > 1 && femcore::par::threads() > 1, |kc| {
            let (mut p_abs, mut p_ext) = (0.0_f64, 0.0_f64);
            let (mut w, mut cu) = (Vec::new(), Vec::new());
            for &ci in chunks[kc] {
                let gid = ned.key_pairs(mesh, ci);
                let nb = ned.nb();
                let (de, _, _) = contrast(&eps[ci], &mu_inv[ci], spec.eps_bg, spec.mu_bg);
                for qp in cell_quadrature(mesh, &geo, ci) {
                    basis3_keyed(ned.order, &qp.lam, &qp.g, &gid, &mut w, &mut cu);
                    let mut es = [c(0.0, 0.0); 3];
                    let mut cs = [c(0.0, 0.0); 3];
                    for s in 0..nb {
                        let coef = x_scat_full[ned.cell_dofs[ci][s]];
                        for k in 0..3 {
                            es[k] += coef * w[s][k];
                            cs[k] += coef * cu[s][k];
                        }
                    }
                    let (e_inc, curl_inc) = prep.field(qp.xq);
                    let e_tot = add_c(es, e_inc);
                    let curl_tot = add_c(cs, curl_inc);
                    // absorption: electric and magnetic loss of the total field
                    let ee = c3_apply(&eps[ci], e_tot);
                    let q_e = e_tot[0].conj() * ee[0] + e_tot[1].conj() * ee[1] + e_tot[2].conj() * ee[2];
                    let mc = c3_apply(&mu_inv[ci], curl_tot);
                    let q_m = curl_tot[0].conj() * mc[0] + curl_tot[1].conj() * mc[1] + curl_tot[2].conj() * mc[2];
                    // Im(curl^H mu_inv curl)/(omega mu0) = omega mu0 mu''|H|^2 (for mu = mu' - i mu'')
                    p_abs += (-0.5 * omega * EPS0 * q_e.im + 0.5 * q_m.im / (omega * MU0)) * qp.wq;
                    if is_contrast[ci] {
                        let dee = c3_apply(&de, e_tot);
                        let q = e_inc[0].conj() * dee[0] + e_inc[1].conj() * dee[1] + e_inc[2].conj() * dee[2];
                        let mut ext = -0.5 * omega * EPS0 * q.im;
                        // magnetic contrast (isotropic mu = 1 / mu_inv[0][0])
                        let mu_c = c(1.0, 0.0) / mu_inv[ci][0][0];
                        let dmu = mu_c - c(spec.mu_bg, 0.0);
                        if dmu.norm() > 1e-12 {
                            // H = i curl E / (omega mu0 mu)
                            let h_tot = scale_c(c(0.0, 1.0 / (omega * MU0)) / mu_c, curl_tot);
                            let h_inc = scale_c(c(0.0, 1.0 / (omega * MU0 * spec.mu_bg)), curl_inc);
                            let qh =
                                h_inc[0].conj() * h_tot[0] + h_inc[1].conj() * h_tot[1] + h_inc[2].conj() * h_tot[2];
                            ext += -0.5 * omega * MU0 * (dmu * qh).im;
                        }
                        p_ext += ext * qp.wq;
                    }
                }
            }
            (p_abs, p_ext)
        });
    let mut p_abs = 0.0;
    let mut p_ext = 0.0;
    for (a, e) in parts {
        p_abs += a;
        p_ext += e;
    }
    (p_abs, p_ext)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curl_fd(prep: &PreparedExcitation, x: [f64; 3], h: f64) -> [C; 3] {
        let mut g = [[c(0.0, 0.0); 3]; 3];
        for a in 0..3 {
            let (mut xp, mut xm) = (x, x);
            xp[a] += h;
            xm[a] -= h;
            let (ep, _) = prep.field(xp);
            let (em, _) = prep.field(xm);
            for k in 0..3 {
                g[a][k] = (ep[k] - em[k]) / (2.0 * h);
            }
        }
        [g[1][2] - g[2][1], g[2][0] - g[0][2], g[0][1] - g[1][0]]
    }

    fn rel(a: [C; 3], b: [C; 3]) -> f64 {
        let d: f64 = (0..3).map(|k| (a[k] - b[k]).norm_sqr()).sum::<f64>().sqrt();
        let n: f64 = (0..3).map(|k| b[k].norm_sqr()).sum::<f64>().sqrt().max(1e-300);
        d / n
    }

    #[test]
    fn analytic_curls_match_finite_differences_and_focus_normalization() {
        let lam = 1.0e-6;
        let omega = 2.0 * PI * C0 / lam;
        let pol = [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)];
        let cases = vec![
            Excitation::PlaneWave { k_hat: [0.0, 0.6, 0.8], pol, amplitude: 2.0, phase: 0.3, origin: [0.0; 3] },
            Excitation::GaussianBeam {
                focus: [0.1e-6, 0.0, 0.0],
                w0: 1.2e-6,
                k_hat: [0.0, 0.0, 1.0],
                pol,
                amplitude: 1.5,
                phase: 0.0,
                model: BeamModel::ComplexSource,
            },
            Excitation::ThinLensBeam {
                focus: [0.0; 3],
                na: 0.6,
                k_hat: [0.0, 0.0, 1.0],
                x_axis: [1.0, 0.0, 0.0],
                pol,
                amplitude: 1.0,
                phase: 0.0,
                filling: 0.0,
                n_theta: 24,
                n_phi: 32,
            },
            Excitation::Dipole {
                position: [0.0; 3],
                moment: [c(0.0, 0.0), c(1e-30, 0.0), c(2e-30, 1e-30)],
                magnetic: false,
            },
            Excitation::Dipole {
                position: [0.0; 3],
                moment: [c(1e-12, 0.0), c(0.0, 0.0), c(0.0, 0.0)],
                magnetic: true,
            },
        ];
        for exc in cases {
            let prep = exc.prepare(omega, c(1.44, 0.0), 1.0).unwrap();
            for x in [[0.31e-6, -0.17e-6, 0.52e-6], [0.9e-6, 0.4e-6, -0.7e-6]] {
                let (_, curl) = prep.field(x);
                let fd = curl_fd(&prep, x, lam * 1e-4);
                assert!(rel(curl, fd) < 1e-5, "{exc:?}: analytic curl vs FD {:.3e}", rel(curl, fd));
            }
        }
        // focus normalization of the complex-source beam and its Gaussian waist
        let w0 = 2.0e-6;
        let gb = Excitation::GaussianBeam {
            focus: [0.0; 3],
            w0,
            k_hat: [0.0, 0.0, 1.0],
            pol,
            amplitude: 1.0,
            phase: 0.0,
            model: BeamModel::ComplexSource,
        };
        let prep = gb.prepare(omega, c(1.0, 0.0), 1.0).unwrap();
        let (e0, _) = prep.field([0.0; 3]);
        assert!((e0[0] - c(1.0, 0.0)).norm() < 1e-12);
        let (e1, _) = prep.field([w0, 0.0, 0.0]);
        // exact beam: |E(w0)| close to 1/e of the paraxial Gaussian (non-paraxial correction ~ 1/(k w0)^2)
        assert!((e1[0].norm() - (-1.0f64).exp()).abs() < 0.02, "|E(w0)| = {}", e1[0].norm());
    }

    #[test]
    fn complex_source_beam_has_no_overflow_for_wide_waists() {
        let lam = 1.0e-6;
        let omega = 2.0 * PI * C0 / lam;
        let gb = Excitation::GaussianBeam {
            focus: [0.0; 3],
            w0: 12.0e-6,
            k_hat: [0.0, 0.0, 1.0],
            pol: [c(0.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)],
            amplitude: 1.0,
            phase: 0.0,
            model: BeamModel::ComplexSource,
        };
        let prep = gb.prepare(omega, c(1.0, 0.0), 1.0).unwrap();
        for x in [[0.0, 0.0, 30e-6], [12e-6, 0.0, 0.0], [12e-6, 5e-6, -40e-6]] {
            let (e, cu) = prep.field(x);
            assert!(e.iter().chain(cu.iter()).all(|z| z.re.is_finite() && z.im.is_finite()));
        }
        let (e, _) = prep.field([0.0, 0.0, 50e-6]);
        // on axis, beyond focus, the amplitude follows 1/sqrt(1 + (z/zR)^2)
        let zr = PI * 12e-6 * 12e-6 / lam;
        let want = 1.0 / (1.0 + (50e-6 / zr).powi(2)).sqrt();
        assert!((e[1].norm() - want).abs() < 5e-3, "on-axis |E| {} vs {}", e[1].norm(), want);
    }
}
