// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Wavelength-dependent optical properties. A material holds dispersion models (constant, Sellmeier, Cauchy, tabulated n and k) per
//! principal axis, evaluated at the study wavelength and temperature: complex index `n - j k` (time convention `e^{j omega t}`, so loss
//! is a negative imaginary part), its derivative with respect to wavelength, the group index and the permittivity tensor.
//! The wavelength unit is declared in the file. A wavelength outside the declared range is an error: nothing extrapolates silently.

use super::reader::{MaterialError, Section};
use super::source::{read_sources, write_sources, Source};
use super::state_t::{clamp_temp_range, AxisClamp, CheckedTempEval};
use super::units::{self, Dim};
use super::M3;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Interp {
    Pchip,
    Linear,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Extrapolate {
    /// outside the table is an error (default)
    Error,
    /// continue with the slope of the end segment, because the file says so
    Linear,
}

/// One dispersion model. Parameters are in the file's wavelength unit (`x` below is the wavelength in that unit).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Dispersion {
    Constant { n: f64, k: f64 },
    /// `n^2 = a + sum_i b_i x^2 / (x^2 - c_i)`, `c_i` in unit^2
    Sellmeier { a: f64, b: Vec<f64>, c: Vec<f64> },
    /// `n = sum_i coefficients[i] x^(-2 i)`
    Cauchy { coefficients: Vec<f64> },
    Tabulated { lambda: Vec<f64>, n: Vec<f64>, k: Vec<f64>, interpolation: Interp, extrapolate: Extrapolate },
    /// Schott glass formula: `n^2 = a0 + a1 x^2 + a2 x^-2 + a3 x^-4 + a4 x^-6 + a5 x^-8`
    Schott { coefficients: Vec<f64> },
    /// Herzberger: `n = A + B L + C L^2 + D x^2 + E x^4 + F x^6`, `L = 1/(x^2 - 0.028 um^2)`
    Herzberger { coefficients: Vec<f64> },
    /// RefractiveIndex.info formula 4 (generalised Sellmeier with power terms, wavelengths in um)
    RiiFormula4 { coefficients: Vec<f64> },
    /// RefractiveIndex.info formula 5 (Cauchy-like general power series `n = c0 + sum c_{1+2i} lam^{c_{2+2i}}`, wavelengths in um)
    RiiFormula5 { coefficients: Vec<f64> },
    /// RefractiveIndex.info formula 6 (gas dispersion `n = 1 + c0 + sum c_{1+2i} / (c_{2+2i} - lam^-2)`, wavelengths in um)
    RiiFormula6 { coefficients: Vec<f64> },
    /// Drude and Lorentz oscillators in the permittivity, `e^{j omega t}` convention:
    /// `eps = eps_inf - f0 Wp^2/(w^2 - j g0 w) + sum_k f_k Wp^2/(w0_k^2 - w^2 + j g_k w)`. Frequencies in rad/s.
    Oscillators { eps_inf: f64, plasma: f64, drude: Option<[f64; 2]>, lorentz: Vec<[f64; 3]>, frequency_unit: String },
}

impl Dispersion {
    pub fn name(&self) -> &'static str {
        match self {
            Dispersion::Constant { .. } => "constant",
            Dispersion::Sellmeier { .. } => "sellmeier",
            Dispersion::Cauchy { .. } => "cauchy",
            Dispersion::Tabulated { .. } => "tabulated",
            Dispersion::Schott { .. } => "schott",
            Dispersion::Herzberger { .. } => "herzberger",
            Dispersion::RiiFormula4 { .. } => "rii_formula_4",
            Dispersion::RiiFormula5 { .. } => "rii_formula_5",
            Dispersion::RiiFormula6 { .. } => "rii_formula_6",
            Dispersion::Oscillators { drude, lorentz, .. } => match (drude.is_some(), lorentz.is_empty()) {
                (true, true) => "drude",
                (false, _) => "lorentz",
                (true, false) => "drude_lorentz",
            },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AxisModel {
    pub dispersion: Dispersion,
    /// where this axis's data come from (if it differs from the `[optical]` table)
    #[serde(default)]
    pub source: Vec<Source>,
    /// thermo-optic coefficient dn/dT in 1/K (0 if none)
    pub dn_dt: f64,
    /// Optional separate extinction table `[wavelengths_in_file_units, k_values]` paired with any dispersion model.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub k_data: Option<(Vec<f64>, Vec<f64>)>,
    /// Optional Einstein oscillator temperature Theta_E [K] saturating dn/dT -> 0 as T -> 0 K.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dn_dt_einstein_theta: Option<f64>,
    /// Optional valid temperature range `[T_min, T_max]` [K] for thermo-optic evaluation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dn_dt_range: Option<[f64; 2]>,
}

impl AxisModel {
    /// Citations attached to this optical axis model.
    pub fn sources(&self) -> &[Source] {
        &self.source
    }

    /// Checks `temperature_k` against `dn_dt_range` (when present) and returns the effective temperature with clamping metadata.
    pub fn eval_dn_dt_temperature_checked(&self, temperature_k: f64) -> CheckedTempEval<f64> {
        match self.dn_dt_range {
            Some([lo, hi]) => {
                let (value, clamp) = clamp_temp_range(temperature_k, lo, hi);
                CheckedTempEval { value, clamp }
            }
            None => CheckedTempEval { value: temperature_k, clamp: None },
        }
    }
}

#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Axes {
    Isotropic(AxisModel),
    /// ordinary for the crystal x and y axes, extraordinary for z
    Uniaxial { ordinary: AxisModel, extraordinary: AxisModel },
    Biaxial([AxisModel; 3]),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Optical {
    /// where the optical data come from
    #[serde(default)]
    pub source: Vec<Source>,
    /// as written in the file (empty if only constant models are used)
    pub wavelength_unit: String,
    /// metres per unit of wavelength
    pub unit_m: f64,
    /// validity range in file units
    pub range: Option<(f64, f64)>,
    /// temperature at which the models hold, K
    pub reference_temperature: f64,
    pub axes: Axes,
}

#[derive(Clone, Debug, PartialEq)]
pub enum OpticsError {
    OutOfRange { lambda_m: f64, lo_m: f64, hi_m: f64, unit: String, unit_m: f64 },
    NotFinite { lambda_m: f64 },
    BadWavelength,
}

impl fmt::Display for OpticsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OpticsError::OutOfRange { lambda_m, lo_m, hi_m, unit, unit_m } => {
                let u = if unit.is_empty() { "m".to_string() } else { unit.clone() };
                let s = if unit.is_empty() { 1.0 } else { *unit_m };
                write!(
                    f,
                    "the wavelength {} {u} is outside the range {} to {} {u} where this material's optical data hold. Choose a wavelength inside the range, or give data that cover it (nothing is extrapolated silently)",
                    lambda_m / s,
                    lo_m / s,
                    hi_m / s
                )
            }
            OpticsError::NotFinite { lambda_m } => write!(f, "the optical model gives no finite value at the wavelength {lambda_m} m"),
            OpticsError::BadWavelength => write!(f, "the wavelength must be a positive, finite number of metres"),
        }
    }
}

impl std::error::Error for OpticsError {}

// ----- evaluation of one model -----

fn pchip_slope(x: &[f64], y: &[f64], i: usize) -> f64 {
    let n = x.len();
    let h = |k: usize| x[k + 1] - x[k];
    let d = |k: usize| (y[k + 1] - y[k]) / h(k);
    if n == 2 {
        return d(0);
    }
    if i == 0 || i == n - 1 {
        let (h0, h1, d0, d1) = if i == 0 { (h(0), h(1), d(0), d(1)) } else { (h(n - 2), h(n - 3), d(n - 2), d(n - 3)) };
        let mut s = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if s * d0 <= 0.0 {
            s = 0.0;
        } else if d0 * d1 <= 0.0 && s.abs() > 3.0 * d0.abs() {
            s = 3.0 * d0;
        }
        return s;
    }
    let (d0, d1) = (d(i - 1), d(i));
    if d0 * d1 <= 0.0 {
        return 0.0;
    }
    let (w1, w2) = (2.0 * h(i) + h(i - 1), h(i) + 2.0 * h(i - 1));
    (w1 + w2) / (w1 / d0 + w2 / d1)
}

/// Value and derivative of a tabulated function at `x`, inside or (with `extend`) outside the table.
fn interpolate(x: &[f64], y: &[f64], interp: Interp, at: f64, extend: bool) -> Option<(f64, f64)> {
    let n = x.len();
    let (lo, hi) = (x[0], x[n - 1]);
    let tol = 1e-12 * (hi - lo).abs().max(hi.abs());
    if at < lo - tol || at > hi + tol {
        if !extend {
            return None;
        }
        let i = if at < lo { 0 } else { n - 1 };
        let slope = match interp {
            Interp::Pchip => pchip_slope(x, y, i),
            Interp::Linear => {
                if i == 0 {
                    (y[1] - y[0]) / (x[1] - x[0])
                } else {
                    (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
                }
            }
        };
        return Some((y[i] + slope * (at - x[i]), slope));
    }
    let at = at.clamp(lo, hi);
    let k = match x.binary_search_by(|v| v.partial_cmp(&at).unwrap()) {
        Ok(k) => k.min(n - 2),
        Err(k) => (k - 1).min(n - 2),
    };
    let h = x[k + 1] - x[k];
    let t = (at - x[k]) / h;
    match interp {
        Interp::Linear => {
            let s = (y[k + 1] - y[k]) / h;
            Some((y[k] + s * (at - x[k]), s))
        }
        Interp::Pchip => {
            let (d0, d1) = (pchip_slope(x, y, k), pchip_slope(x, y, k + 1));
            let (t2, t3) = (t * t, t * t * t);
            let v = (2.0 * t3 - 3.0 * t2 + 1.0) * y[k] + (t3 - 2.0 * t2 + t) * h * d0 + (-2.0 * t3 + 3.0 * t2) * y[k + 1] + (t3 - t2) * h * d1;
            let dv = ((6.0 * t2 - 6.0 * t) * y[k] + (3.0 * t2 - 4.0 * t + 1.0) * h * d0 + (-6.0 * t2 + 6.0 * t) * y[k + 1] + (3.0 * t2 - 2.0 * t) * h * d1) / h;
            Some((v, dv))
        }
    }
}

impl Dispersion {
    /// `(n, k, dn/dx, dk/dx)` at the wavelength `x` (file units); `None` if a table does not cover `x`.
    fn eval(&self, x: f64, unit_m: f64) -> Option<(f64, f64, f64, f64)> {
        match self {
            Dispersion::Constant { n, k } => Some((*n, *k, 0.0, 0.0)),
            Dispersion::Sellmeier { a, b, c } => {
                let x2 = x * x;
                let mut n2 = *a;
                let mut dn2 = 0.0;
                for (bi, ci) in b.iter().zip(c) {
                    let den = x2 - ci;
                    n2 += bi * x2 / den;
                    dn2 += -2.0 * bi * ci * x / (den * den);
                }
                if n2 <= 0.0 || !n2.is_finite() {
                    return None;
                }
                let n = n2.sqrt();
                Some((n, 0.0, dn2 / (2.0 * n), 0.0))
            }
            Dispersion::Cauchy { coefficients } => {
                let (mut n, mut dn) = (0.0, 0.0);
                for (i, c) in coefficients.iter().enumerate() {
                    let p = 2.0 * i as f64;
                    n += c * x.powf(-p);
                    dn += -p * c * x.powf(-p - 1.0);
                }
                Some((n, 0.0, dn, 0.0))
            }
            Dispersion::Tabulated { lambda, n, k, interpolation, extrapolate } => {
                let ext = *extrapolate == Extrapolate::Linear;
                let (nv, dn) = interpolate(lambda, n, *interpolation, x, ext)?;
                let (kv, dk) = interpolate(lambda, k, *interpolation, x, ext)?;
                Some((nv, kv, dn, dk))
            }
            Dispersion::Schott { coefficients: a } => {
                let n2 = a[0] + a[1] * x.powi(2) + a[2] * x.powi(-2) + a[3] * x.powi(-4) + a[4] * x.powi(-6) + a[5] * x.powi(-8);
                let dn2 = 2.0 * a[1] * x - 2.0 * a[2] * x.powi(-3) - 4.0 * a[3] * x.powi(-5) - 6.0 * a[4] * x.powi(-7) - 8.0 * a[5] * x.powi(-9);
                if n2 <= 0.0 || !n2.is_finite() {
                    return None;
                }
                let n = n2.sqrt();
                Some((n, 0.0, dn2 / (2.0 * n), 0.0))
            }
            Dispersion::Herzberger { coefficients: c } => {
                // the constant 0.028 is in um^2: express it in the square of the file's wavelength unit
                let l0 = 0.028e-12 / (unit_m * unit_m);
                let d = x * x - l0;
                if d == 0.0 {
                    return None;
                }
                let n = c[0] + c[1] / d + c[2] / (d * d) + c[3] * x.powi(2) + c[4] * x.powi(4) + c[5] * x.powi(6);
                let dn = -2.0 * x * c[1] / (d * d) - 4.0 * x * c[2] / (d * d * d) + 2.0 * c[3] * x + 4.0 * c[4] * x.powi(3) + 6.0 * c[5] * x.powi(5);
                if !n.is_finite() {
                    return None;
                }
                Some((n, 0.0, dn, 0.0))
            }
            Dispersion::RiiFormula4 { coefficients: c } => {
                if c.is_empty() {
                    return None;
                }
                let s_um = unit_m / 1e-6;
                let xu = x * s_um;
                let x2 = xu * xu;
                let mut n2 = c[0];
                let mut dn2 = 0.0;
                let mut i = 1;
                for _ in 0..2 {
                    if i + 3 < c.len() {
                        if c[i] != 0.0 {
                            let pole = c[i + 2].powf(c[i + 3]);
                            let den = x2 - pole;
                            if den.abs() < 1e-18 || !den.is_finite() {
                                return None;
                            }
                            let p = c[i + 1];
                            let xp = xu.powf(p);
                            let num = c[i] * xp;
                            let dnum = if p == 0.0 { 0.0 } else { c[i] * p * xu.powf(p - 1.0) };
                            n2 += num / den;
                            dn2 += (dnum * den - num * (2.0 * xu)) / (den * den);
                        }
                        i += 4;
                    }
                }
                while i + 1 < c.len() {
                    let a = c[i];
                    let p = c[i + 1];
                    n2 += a * xu.powf(p);
                    if p != 0.0 {
                        dn2 += a * p * xu.powf(p - 1.0);
                    }
                    i += 2;
                }
                if n2 <= 0.0 || !n2.is_finite() {
                    return None;
                }
                let n = n2.sqrt();
                Some((n, 0.0, (dn2 / (2.0 * n)) * s_um, 0.0))
            }
            Dispersion::RiiFormula5 { coefficients: c } => {
                if c.is_empty() {
                    return None;
                }
                let s_um = unit_m / 1e-6;
                let xu = x * s_um;
                let mut n = c[0];
                let mut dn = 0.0;
                for p in c[1..].as_chunks::<2>().0 {
                    n += p[0] * xu.powf(p[1]);
                    if p[1] != 0.0 {
                        dn += p[0] * p[1] * xu.powf(p[1] - 1.0);
                    }
                }
                if n <= 0.0 || !n.is_finite() {
                    return None;
                }
                Some((n, 0.0, dn * s_um, 0.0))
            }
            Dispersion::RiiFormula6 { coefficients: c } => {
                if c.is_empty() {
                    return None;
                }
                let s_um = unit_m / 1e-6;
                let xu = x * s_um;
                let inv = 1.0 / (xu * xu);
                let inv3 = 1.0 / (xu * xu * xu);
                let mut n = 1.0 + c[0];
                let mut dn = 0.0;
                for p in c[1..].as_chunks::<2>().0 {
                    let den = p[1] - inv;
                    if den.abs() < 1e-18 || !den.is_finite() {
                        return None;
                    }
                    n += p[0] / den;
                    dn += -2.0 * p[0] * inv3 / (den * den);
                }
                if n <= 0.0 || !n.is_finite() {
                    return None;
                }
                Some((n, 0.0, dn * s_um, 0.0))
            }
            Dispersion::Oscillators { eps_inf, plasma, drude, lorentz, .. } => {
                use num_complex::Complex64 as C;
                let lam = x * unit_m;
                let w = 2.0 * std::f64::consts::PI * 299_792_458.0 / lam;
                let j = C::new(0.0, 1.0);
                let wp2 = plasma * plasma;
                let mut eps = C::new(*eps_inf, 0.0);
                let mut deps = C::new(0.0, 0.0); // d eps / d w
                if let Some([f0, g0]) = drude {
                    let q = C::new(w * w, 0.0) - j * (g0 * w);
                    eps -= (f0 * wp2) / q;
                    deps += (f0 * wp2) * (C::new(2.0 * w, 0.0) - j * *g0) / (q * q);
                }
                for [f, w0, g] in lorentz {
                    let dd = C::new(w0 * w0 - w * w, 0.0) + j * (g * w);
                    eps += (f * wp2) / dd;
                    deps += (f * wp2) * (C::new(2.0 * w, 0.0) - j * *g) / (dd * dd);
                }
                let nc = eps.sqrt();
                if !(nc.re.is_finite() && nc.im.is_finite()) || nc.norm() == 0.0 {
                    return None;
                }
                let dn_dlam = deps / (nc * 2.0) * (-w / lam);
                // n - j k = sqrt(eps), so k = -Im; derivatives are per file unit
                Some((nc.re, -nc.im, dn_dlam.re * unit_m, -dn_dlam.im * unit_m))
            }
        }
    }

    /// Smallest and largest wavelength a table covers (file units); `None` for formulas.
    fn table_extent(&self) -> Option<(f64, f64)> {
        match self {
            Dispersion::Tabulated { lambda, extrapolate: Extrapolate::Error, .. } => Some((lambda[0], lambda[lambda.len() - 1])),
            _ => None,
        }
    }
}

impl Axes {
    pub fn model(&self, axis: usize) -> &AxisModel {
        match self {
            Axes::Isotropic(m) => m,
            Axes::Uniaxial { ordinary, extraordinary } => {
                if axis == 2 {
                    extraordinary
                } else {
                    ordinary
                }
            }
            Axes::Biaxial(m) => &m[axis.min(2)],
        }
    }
    pub fn models(&self) -> Vec<&AxisModel> {
        match self {
            Axes::Isotropic(m) => vec![m],
            Axes::Uniaxial { ordinary, extraordinary } => vec![ordinary, extraordinary],
            Axes::Biaxial(m) => m.iter().collect(),
        }
    }
}

impl Optical {
    pub fn is_anisotropic(&self) -> bool {
        !matches!(self.axes, Axes::Isotropic(_))
    }

    /// The validity range in metres.
    pub fn range_m(&self) -> Option<(f64, f64)> {
        self.range.map(|(lo, hi)| (lo * self.unit_m, hi * self.unit_m))
    }

    /// `(n, k, dn/dlambda, dk/dlambda)` for crystal axis 0, 1 or 2 (derivatives per metre), returning any `dn_dT_range` clamping metadata.
    pub fn nk_checked(
        &self,
        axis: usize,
        lambda_m: f64,
        temp_k: Option<f64>,
    ) -> Result<CheckedTempEval<(f64, f64, f64, f64)>, OpticsError> {
        if !(lambda_m.is_finite() && lambda_m > 0.0) {
            return Err(OpticsError::BadWavelength);
        }
        let x = lambda_m / self.unit_m;
        let m = self.axes.model(axis);
        let out_of = |lo: f64, hi: f64| OpticsError::OutOfRange { lambda_m, lo_m: lo * self.unit_m, hi_m: hi * self.unit_m, unit: self.wavelength_unit.clone(), unit_m: self.unit_m };
        let tol = 1e-12 * x.abs();
        if let Some((lo, hi)) = self.range {
            if x < lo - tol || x > hi + tol {
                return Err(out_of(lo, hi));
            }
        }
        if let Some((lo, hi)) = m.dispersion.table_extent() {
            if x < lo - tol || x > hi + tol {
                return Err(out_of(lo, hi));
            }
        }
        let (mut n, mut k, dn, mut dk) = m.dispersion.eval(x, self.unit_m).ok_or(OpticsError::NotFinite { lambda_m })?;
        if let Some((kl, kv)) = &m.k_data {
            if let Some((k_val, dk_val)) = interpolate(kl, kv, Interp::Pchip, x, false) {
                k = k_val.max(0.0);
                dk = dk_val;
            } else if x > kl[kl.len() - 1] && kv[kv.len() - 1] <= 1e-4 {
                // Transparent long-wavelength tail past interband k(lambda) table
                k = 0.0;
                dk = 0.0;
            } else {
                return Err(out_of(kl[0], kl[kl.len() - 1]));
            }
        }
        let mut clamp = None;
        if let Some(t_raw) = temp_k {
            let tc = m.eval_dn_dt_temperature_checked(t_raw);
            clamp = tc.clamp;
            let t = tc.value;
            if let Some(theta) = m.dn_dt_einstein_theta {
                let g_t = super::state_t::einstein_thermal_coordinate(t, theta);
                let g_ref = super::state_t::einstein_thermal_coordinate(self.reference_temperature, theta);
                n += m.dn_dt * (g_t - g_ref);
            } else {
                n += m.dn_dt * (t - self.reference_temperature);
            }
        }
        if !(n.is_finite() && k.is_finite()) {
            return Err(OpticsError::NotFinite { lambda_m });
        }
        Ok(CheckedTempEval { value: (n, k, dn / self.unit_m, dk / self.unit_m), clamp })
    }

    /// `(n, k, dn/dlambda, dk/dlambda)` for crystal axis 0, 1 or 2 (derivatives per metre).
    pub fn nk(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<(f64, f64, f64, f64), OpticsError> {
        self.nk_checked(axis, lambda_m, temp_k).map(|r| r.value)
    }

    /// Complex index `n - j k` together with any `dn_dT_range` clamping metadata.
    pub fn index_checked(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<CheckedTempEval<Complex64>, OpticsError> {
        let checked = self.nk_checked(axis, lambda_m, temp_k)?;
        Ok(CheckedTempEval {
            value: Complex64::new(checked.value.0, -checked.value.1),
            clamp: checked.clamp,
        })
    }

    /// Complex index `n - j k`.
    pub fn index(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<Complex64, OpticsError> {
        self.index_checked(axis, lambda_m, temp_k).map(|r| r.value)
    }

    /// `d(n - j k)/d lambda` per metre.
    pub fn d_index_d_lambda(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<Complex64, OpticsError> {
        let (_, _, dn, dk) = self.nk(axis, lambda_m, temp_k)?;
        Ok(Complex64::new(dn, -dk))
    }

    /// Group index `n - lambda dn/dlambda` (real part).
    pub fn group_index(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<f64, OpticsError> {
        let (n, _, dn, _) = self.nk(axis, lambda_m, temp_k)?;
        Ok(n - lambda_m * dn)
    }

    /// Relative permittivity along the three crystal axes, `(n - j k)^2`, returning any per-axis clamping events.
    pub fn permittivity_diag_checked(
        &self,
        lambda_m: f64,
        temp_k: Option<f64>,
    ) -> Result<([Complex64; 3], Vec<AxisClamp>), OpticsError> {
        let mut e = [Complex64::new(0.0, 0.0); 3];
        let mut clamps = Vec::new();
        for (a, v) in e.iter_mut().enumerate() {
            let checked = self.index_checked(a, lambda_m, temp_k)?;
            if let Some(info) = checked.clamp {
                clamps.push((a, info));
            }
            *v = checked.value * checked.value;
        }
        Ok((e, clamps))
    }

    /// Relative permittivity along the three crystal axes, `(n - j k)^2`.
    pub fn permittivity_diag(&self, lambda_m: f64, temp_k: Option<f64>) -> Result<[Complex64; 3], OpticsError> {
        self.permittivity_diag_checked(lambda_m, temp_k).map(|(e, _)| e)
    }

    /// Relative permittivity tensor in the global frame for a rotated crystal, returning any per-axis `dn_dT_range` clamping events.
    pub fn permittivity_global_checked(
        &self,
        lambda_m: f64,
        temp_k: Option<f64>,
        euler_deg: [f64; 3],
    ) -> Result<([[Complex64; 3]; 3], Vec<AxisClamp>), OpticsError> {
        let (d, clamps) = self.permittivity_diag_checked(lambda_m, temp_k)?;
        let r: M3 = super::rotation(euler_deg);
        let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for (k, dk) in d.iter().enumerate() {
                    out[i][j] += *dk * (r[i][k] * r[j][k]);
                }
            }
        }
        Ok((out, clamps))
    }

    /// Relative permittivity tensor in the global frame for a crystal rotated by the Bunge angles (degrees).
    pub fn permittivity_global(&self, lambda_m: f64, temp_k: Option<f64>, euler_deg: [f64; 3]) -> Result<[[Complex64; 3]; 3], OpticsError> {
        self.permittivity_global_checked(lambda_m, temp_k, euler_deg).map(|(m, _)| m)
    }

    // ----- reading and writing -----

    pub fn read(s: &mut Section) -> Result<Optical, MaterialError> {
        let unit_text = s.string("wavelength_unit")?;
        let unit_m = match &unit_text {
            Some(u) => units::to_si(1.0, u, Dim::Length).map_err(|m| s.err("wavelength_unit", m))?,
            None => 1.0,
        };
        let range = match s.value("range") {
            None => None,
            Some(toml::Value::Array(a)) if a.len() == 2 => {
                let mut r = [0.0; 2];
                for (i, v) in a.iter().enumerate() {
                    r[i] = match v {
                        toml::Value::String(q) => {
                            let (x, u) = units::split_quantity(q).map_err(|m| s.err("range", m))?;
                            if u.is_empty() {
                                x
                            } else {
                                units::to_si(x, &u, Dim::Length).map_err(|m| s.err("range", m))? / unit_m
                            }
                        }
                        other => s.number_of("range", other)?,
                    };
                }
                if unit_text.is_none() && a.iter().any(|v| !v.is_str()) {
                    return Err(s.err("wavelength_unit", "is required: the range and the model coefficients are in this unit (for example wavelength_unit = \"um\"); units are never assumed"));
                }
                if !(r[0] > 0.0 && r[1] > r[0]) {
                    return Err(s.err("range", "must be two wavelengths, the smaller first, both above zero, for example [0.21, 6.7]"));
                }
                Some((r[0], r[1]))
            }
            Some(_) => return Err(s.err("range", "must be a list of two wavelengths, for example [0.21, 6.7]")),
        };
        let reference_temperature = s.quantity("reference_temperature", Dim::Temperature, None)?.unwrap_or(293.15);
        let source = read_sources(s, "source")?;

        let sub_names = ["ordinary", "extraordinary", "x", "y", "z"];
        let has_sub: Vec<&str> = sub_names.iter().copied().filter(|k| s.has(k)).collect();
        let ctx = Ctx { unit_m, unit_given: unit_text.is_some(), range };
        let axes = if has_sub.is_empty() {
            Axes::Isotropic(read_model(s, &ctx)?)
        } else {
            if s.has("model") {
                return Err(s.err("model", "cannot be combined with the per-axis tables [optical.ordinary]/[optical.extraordinary] or [optical.x]/[.y]/[.z]; put the model inside each table"));
            }
            let uni = has_sub.iter().all(|k| *k == "ordinary" || *k == "extraordinary");
            let bi = has_sub.iter().all(|k| *k == "x" || *k == "y" || *k == "z");
            if uni {
                if has_sub.len() != 2 {
                    return Err(s.err("", "a uniaxial crystal needs both [optical.ordinary] and [optical.extraordinary]"));
                }
                let mut o = s.sub("ordinary")?.unwrap();
                let mut e = s.sub("extraordinary")?.unwrap();
                Axes::Uniaxial { ordinary: read_model_finish(&mut o, &ctx)?, extraordinary: read_model_finish(&mut e, &ctx)? }
            } else if bi {
                if has_sub.len() != 3 {
                    return Err(s.err("", "a biaxial crystal needs all of [optical.x], [optical.y] and [optical.z]"));
                }
                let mut m = Vec::new();
                for k in ["x", "y", "z"] {
                    let mut t = s.sub(k)?.unwrap();
                    m.push(read_model_finish(&mut t, &ctx)?);
                }
                Axes::Biaxial([m[0].clone(), m[1].clone(), m[2].clone()])
            } else {
                return Err(s.err("", "mixes ordinary/extraordinary with x/y/z tables; use one of the two forms"));
            }
        };
        let opt = Optical { source, wavelength_unit: unit_text.unwrap_or_default(), unit_m, range, reference_temperature, axes };
        opt.validate().map_err(|(path, m)| s.err(&path, m))?;
        Ok(opt)
    }

    /// Checks that every model is usable over the declared range.
    pub fn validate(&self) -> Result<(), (String, String)> {
        for (idx, m) in self.axes.models().into_iter().enumerate() {
            let tag = match &self.axes {
                Axes::Isotropic(_) => String::new(),
                Axes::Uniaxial { .. } => ["ordinary", "extraordinary"][idx].to_string(),
                Axes::Biaxial(_) => ["x", "y", "z"][idx].to_string(),
            };
            let at = |m: String| (tag.clone(), m);
            let Some((lo, hi)) = self.range.or_else(|| m.dispersion.table_extent()) else {
                if matches!(m.dispersion, Dispersion::Constant { .. }) {
                    continue;
                }
                return Err(at(format!(
                    "the {} model needs a validity range, for example range = [0.21, 6.7] in the wavelength unit: published coefficients hold only there and nothing is extrapolated silently",
                    m.dispersion.name()
                )));
            };
            match &m.dispersion {
                Dispersion::Constant { n, k } => {
                    if *n <= 0.0 || *k < 0.0 {
                        return Err(at("n must be positive and k must not be negative".into()));
                    }
                }
                Dispersion::Sellmeier { b, c, .. } => {
                    for ci in c {
                        let pole = ci.sqrt();
                        if pole >= lo && pole <= hi {
                            return Err(at(format!("the Sellmeier resonance at {pole} (unit of the wavelength) lies inside the range {lo} to {hi}, where n would diverge; check C or the range")));
                        }
                    }
                    let _ = b;
                }
                Dispersion::Cauchy { .. }
                | Dispersion::Schott { .. }
                | Dispersion::Herzberger { .. }
                | Dispersion::RiiFormula4 { .. }
                | Dispersion::RiiFormula5 { .. }
                | Dispersion::RiiFormula6 { .. }
                | Dispersion::Oscillators { .. } => {}
                Dispersion::Tabulated { lambda, n, k, .. } => {
                    if let Some((tlo, thi)) = self.range {
                        if m.dispersion.table_extent().is_some() && (tlo < lambda[0] * (1.0 - 1e-12) || thi > lambda[lambda.len() - 1] * (1.0 + 1e-12)) {
                            return Err(at(format!(
                                "the range {tlo} to {thi} is wider than the table ({} to {}); shorten the range, extend the table, or write extrapolate = \"linear\" to allow it on purpose",
                                lambda[0],
                                lambda[lambda.len() - 1]
                            )));
                        }
                    }
                    for (i, (nv, kv)) in n.iter().zip(k).enumerate() {
                        if *nv <= 0.0 {
                            return Err(at(format!("table row {}: n must be positive (found {nv})", i + 1)));
                        }
                        if *kv < 0.0 {
                            return Err(at(format!("table row {}: k must not be negative (found {kv}); loss is written as a positive k", i + 1)));
                        }
                    }
                }
            }
            // sample the model over the range
            if !matches!(m.dispersion, Dispersion::Constant { .. } | Dispersion::Tabulated { .. }) {
                for i in 0..=128 {
                    let x = lo + (hi - lo) * i as f64 / 128.0;
                    match m.dispersion.eval(x, self.unit_m) {
                        Some((n, k, ..)) if n.is_finite() && k.is_finite() && (n > 0.0 || (n >= 0.0 && matches!(m.dispersion, Dispersion::Oscillators { .. }))) => {}
                        _ => return Err(at(format!("the {} model gives no positive, finite index at the wavelength {x} (inside the declared range {lo} to {hi}); check the coefficients and the unit", m.dispersion.name()))),
                    }
                }
            }
        }
        Ok(())
    }

    pub fn write(&self) -> toml::Table {
        let mut t = toml::Table::new();
        if let Some(v) = write_sources(&self.source) {
            t.insert("source".into(), v);
        }
        if !self.wavelength_unit.is_empty() {
            t.insert("wavelength_unit".into(), self.wavelength_unit.clone().into());
        }
        if let Some((lo, hi)) = self.range {
            t.insert("range".into(), toml::Value::Array(vec![toml::Value::Float(lo), toml::Value::Float(hi)]));
        }
        t.insert("reference_temperature".into(), format!("{} K", self.reference_temperature).into());
        match &self.axes {
            Axes::Isotropic(m) => write_model(&mut t, m),
            Axes::Uniaxial { ordinary, extraordinary } => {
                for (k, m) in [("ordinary", ordinary), ("extraordinary", extraordinary)] {
                    let mut s = toml::Table::new();
                    write_model(&mut s, m);
                    t.insert(k.into(), toml::Value::Table(s));
                }
            }
            Axes::Biaxial(ms) => {
                for (k, m) in ["x", "y", "z"].into_iter().zip(ms) {
                    let mut s = toml::Table::new();
                    write_model(&mut s, m);
                    t.insert(k.into(), toml::Value::Table(s));
                }
            }
        }
        t
    }
}

struct Ctx {
    unit_m: f64,
    unit_given: bool,
    range: Option<(f64, f64)>,
}

fn read_model_finish(s: &mut Section, ctx: &Ctx) -> Result<AxisModel, MaterialError> {
    let mut m = read_model(s, ctx)?;
    m.source = read_sources(s, "source")?;
    s.finish(MODEL_KEYS)?;
    Ok(m)
}

pub const MODEL_KEYS: &[&str] = &[
    "model",
    "n",
    "k",
    "eps",
    "A",
    "B",
    "C",
    "lambda0",
    "coefficients",
    "data",
    "k_data",
    "interpolation",
    "extrapolate",
    "dn_dT",
    "dn_dT_unit",
    "dn_dT_einstein_theta",
    "dn_dT_einstein_theta_unit",
    "dn_dT_range",
    "dn_dT_range_unit",
    "source",
    "frequency_unit",
    "eps_inf",
    "plasma_frequency",
    "damping",
    "drude_weight",
    "oscillators",
];

fn floats(s: &mut Section, key: &str) -> Result<Option<Vec<f64>>, MaterialError> {
    s.numbers(key)
}

type KDataTable = (Vec<f64>, Vec<f64>);

fn read_k_data(s: &mut Section, need_unit: bool) -> Result<Option<KDataTable>, MaterialError> {
    let Some(v) = s.value("k_data") else {
        return Ok(None);
    };
    if !need_unit {
        return Err(s.err("wavelength_unit", "is required in [optical] when using k_data"));
    }
    let toml::Value::Array(rows) = v else {
        return Err(s.err("k_data", "must be [[wavelength, k], ...]"));
    };
    let mut kl = Vec::with_capacity(rows.len());
    let mut kv = Vec::with_capacity(rows.len());
    for (i, row) in rows.iter().enumerate() {
        let toml::Value::Array(r) = row else {
            return Err(s.err("k_data", format!("row {} must be [wavelength, k]", i + 1)));
        };
        if r.len() != 2 {
            return Err(s.err("k_data", format!("row {} has {} entries; expected [wavelength, k]", i + 1, r.len())));
        }
        let wl = s.number_of("k_data", &r[0])?;
        let k = s.number_of("k_data", &r[1])?;
        if wl <= 0.0 || k < 0.0 {
            return Err(s.err("k_data", format!("row {}: wavelength must be > 0 and k must be >= 0", i + 1)));
        }
        if kl.last().is_some_and(|&prev: &f64| wl <= prev) {
            return Err(s.err("k_data", format!("row {}: wavelengths must strictly increase", i + 1)));
        }
        kl.push(wl);
        kv.push(k);
    }
    if kl.len() < 2 {
        return Err(s.err("k_data", "k_data requires at least two [wavelength, k] rows"));
    }
    Ok(Some((kl, kv)))
}

fn read_dn_dt_range(s: &mut Section) -> Result<Option<[f64; 2]>, MaterialError> {
    let Some(v) = s.value("dn_dT_range") else {
        return Ok(None);
    };
    let u = s.string("dn_dT_range_unit")?.unwrap_or_else(|| "K".to_string());
    let toml::Value::Array(a) = v else {
        return Err(s.err("dn_dT_range", "must be [T_min, T_max]"));
    };
    if a.len() != 2 {
        return Err(s.err("dn_dT_range", "dn_dT_range must have 2 elements [T_min, T_max]"));
    }
    let lo_raw = s.number_of("dn_dT_range", &a[0])?;
    let hi_raw = s.number_of("dn_dT_range", &a[1])?;
    let lo = units::to_si(lo_raw, &u, Dim::Temperature).map_err(|m| s.err("dn_dT_range", m))?;
    let hi = units::to_si(hi_raw, &u, Dim::Temperature).map_err(|m| s.err("dn_dT_range", m))?;
    if lo < 0.0 || hi <= lo {
        return Err(s.err("dn_dT_range", "dn_dT_range must satisfy 0 <= T_min < T_max"));
    }
    Ok(Some([lo, hi]))
}

fn read_model(s: &mut Section, ctx: &Ctx) -> Result<AxisModel, MaterialError> {
    let kind = s.string("model")?.ok_or_else(|| s.err("model", "is missing; write model = \"constant\", \"sellmeier\", \"cauchy\" or \"tabulated\""))?;
    let dn_dt = s.quantity("dn_dT", Dim::PerKelvin, None)?.unwrap_or(0.0);
    let dn_dt_einstein_theta = s.quantity("dn_dT_einstein_theta", Dim::Temperature, Some("K"))?;
    if let Some(th) = dn_dt_einstein_theta {
        if th <= 0.0 {
            return Err(s.err("dn_dT_einstein_theta", "Einstein temperature Theta_E must be positive"));
        }
    }
    let dn_dt_range = read_dn_dt_range(s)?;
    let k_data = read_k_data(s, ctx.unit_given)?;
    let need_unit = |s: &Section| -> Result<(), MaterialError> {
        if ctx.unit_given {
            Ok(())
        } else {
            Err(s.err("wavelength_unit", "is required in [optical] for this model: its coefficients are in this unit (for example wavelength_unit = \"um\")"))
        }
    };
    let dispersion = match kind.to_ascii_lowercase().as_str() {
        "constant" => {
            let n = match (s.number("n")?, s.number("eps")?) {
                (Some(n), None) => n,
                (None, Some(e)) => {
                    if e <= 0.0 {
                        return Err(s.err("eps", "the relative permittivity of an optical material must be positive; absorbing media use k"));
                    }
                    e.sqrt()
                }
                (None, None) => return Err(s.err("", "a constant model needs n (refractive index) or eps (relative permittivity)")),
                _ => return Err(s.err("", "give n or eps, not both")),
            };
            let k = s.number("k")?.unwrap_or(0.0);
            Dispersion::Constant { n, k }
        }
        "sellmeier" => {
            need_unit(s)?;
            let a = s.number("A")?.unwrap_or(1.0);
            let b = floats(s, "B")?.ok_or_else(|| s.err("B", "the Sellmeier model needs B = [B1, B2, ...] and C = [C1, C2, ...] (C in wavelength_unit squared) or lambda0 = [...] (resonance wavelengths)"))?;
            let c = match (floats(s, "C")?, floats(s, "lambda0")?) {
                (Some(c), None) => c,
                (None, Some(l)) => l.iter().map(|x| x * x).collect(),
                (None, None) => return Err(s.err("C", "is missing; give C (squared resonance wavelengths) or lambda0 (resonance wavelengths), in the wavelength unit")),
                _ => return Err(s.err("C", "give C or lambda0, not both")),
            };
            if b.len() != c.len() || b.is_empty() {
                return Err(s.err("B", format!("B has {} entries but C has {}; there must be one of each per term, and at least one", b.len(), c.len())));
            }
            if c.iter().any(|x| *x <= 0.0) {
                return Err(s.err("C", "every C must be positive (it is the square of a resonance wavelength)"));
            }
            Dispersion::Sellmeier { a, b, c }
        }
        "cauchy" => {
            need_unit(s)?;
            let coefficients = floats(s, "coefficients")?.ok_or_else(|| s.err("coefficients", "the Cauchy model needs coefficients = [A, B, C, ...] for n = A + B/x^2 + C/x^4 + ..."))?;
            if coefficients.is_empty() {
                return Err(s.err("coefficients", "must have at least one entry"));
            }
            Dispersion::Cauchy { coefficients }
        }
        "tabulated" => {
            need_unit(s)?;
            let rows = match s.value("data") {
                Some(toml::Value::Array(r)) => r,
                _ => return Err(s.err("data", "a tabulated model needs data = [[wavelength, n], ...] or [[wavelength, n, k], ...]")),
            };
            let (mut lambda, mut n, mut k) = (Vec::new(), Vec::new(), Vec::new());
            let mut width = 0;
            for (i, row) in rows.iter().enumerate() {
                let toml::Value::Array(r) = row else { return Err(s.err("data", format!("row {} must be a list [wavelength, n] or [wavelength, n, k]", i + 1))) };
                if r.len() != 2 && r.len() != 3 {
                    return Err(s.err("data", format!("row {} has {} entries; expected [wavelength, n] or [wavelength, n, k]", i + 1, r.len())));
                }
                if width != 0 && width != r.len() {
                    return Err(s.err("data", format!("row {} has {} entries but the rows before have {}; all rows need the same columns", i + 1, r.len(), width)));
                }
                width = r.len();
                let v: Vec<f64> = r.iter().map(|x| s.number_of("data", x)).collect::<Result<_, _>>()?;
                lambda.push(v[0]);
                n.push(v[1]);
                k.push(if v.len() == 3 { v[2] } else { 0.0 });
            }
            if lambda.len() < 2 {
                return Err(s.err("data", "needs at least two rows"));
            }
            for i in 1..lambda.len() {
                if lambda[i] <= lambda[i - 1] {
                    return Err(s.err("data", format!("the wavelengths must increase: row {} has {} after {}", i + 1, lambda[i], lambda[i - 1])));
                }
            }
            if lambda[0] <= 0.0 {
                return Err(s.err("data", "wavelengths must be positive"));
            }
            let interpolation = match s.string("interpolation")?.as_deref() {
                None | Some("pchip") => Interp::Pchip,
                Some("linear") => Interp::Linear,
                Some(o) => return Err(s.err("interpolation", format!("'{o}' is not known; use \"pchip\" (default, monotone cubic) or \"linear\""))),
            };
            let extrapolate = match s.string("extrapolate")?.as_deref() {
                None | Some("error") => Extrapolate::Error,
                Some("linear") => Extrapolate::Linear,
                Some(o) => return Err(s.err("extrapolate", format!("'{o}' is not known; use \"error\" (default) or \"linear\""))),
            };
            if extrapolate == Extrapolate::Error {
                if let Some((lo, hi)) = ctx.range {
                    if lo < lambda[0] * (1.0 - 1e-12) || hi > lambda[lambda.len() - 1] * (1.0 + 1e-12) {
                        return Err(s.err("range", format!("the range {lo} to {hi} is wider than the table ({} to {}); shorten the range, extend the table, or write extrapolate = \"linear\" on purpose", lambda[0], lambda[lambda.len() - 1])));
                    }
                }
            }
            let _ = ctx.unit_m;
            Dispersion::Tabulated { lambda, n, k, interpolation, extrapolate }
        }
        "schott" | "herzberger" | "rii_formula_4" | "rii_formula_5" | "rii_formula_6" => {
            need_unit(s)?;
            let name = kind.to_ascii_lowercase();
            let coefficients = floats(s, "coefficients")?.ok_or_else(|| {
                match name.as_str() {
                    "schott" => s.err("coefficients", "the Schott model needs coefficients = [a0, a1, a2, a3, a4, a5] for n^2 = a0 + a1 x^2 + a2 x^-2 + a3 x^-4 + a4 x^-6 + a5 x^-8"),
                    "herzberger" => s.err("coefficients", "the Herzberger model needs coefficients = [A, B, C, D, E, F] for n = A + B L + C L^2 + D x^2 + E x^4 + F x^6, L = 1/(x^2 - 0.028 um^2)"),
                    _ => s.err("coefficients", format!("the {name} model requires a coefficients = [...] list")),
                }
            })?;
            match name.as_str() {
                "schott" => {
                    if coefficients.len() != 6 {
                        return Err(s.err("coefficients", format!("the schott model needs exactly 6 coefficients (found {}); write 0 for the ones you do not have", coefficients.len())));
                    }
                    Dispersion::Schott { coefficients }
                }
                "herzberger" => {
                    if coefficients.len() != 6 {
                        return Err(s.err("coefficients", format!("the herzberger model needs exactly 6 coefficients (found {}); write 0 for the ones you do not have", coefficients.len())));
                    }
                    Dispersion::Herzberger { coefficients }
                }
                "rii_formula_4" => {
                    if coefficients.is_empty() {
                        return Err(s.err("coefficients", "rii_formula_4 needs at least 1 coefficient"));
                    }
                    Dispersion::RiiFormula4 { coefficients }
                }
                "rii_formula_5" => {
                    if coefficients.is_empty() || coefficients.len() % 2 != 1 {
                        return Err(s.err("coefficients", "rii_formula_5 needs 1 + 2*m coefficients [c0, c1, c2, ...]"));
                    }
                    Dispersion::RiiFormula5 { coefficients }
                }
                _ => {
                    if coefficients.is_empty() || coefficients.len() % 2 != 1 {
                        return Err(s.err("coefficients", "rii_formula_6 needs 1 + 2*m coefficients [c0, c1, c2, ...]"));
                    }
                    Dispersion::RiiFormula6 { coefficients }
                }
            }
        }
        "drude" | "lorentz" | "drude_lorentz" => {
            let name = kind.to_ascii_lowercase();
            let fu = s
                .string("frequency_unit")?
                .ok_or_else(|| s.err("frequency_unit", "is required for this model: the frequencies are in this unit (for example frequency_unit = \"eV\", \"THz\" or \"rad/s\")"))?;
            let scale = units::to_si(1.0, &fu, Dim::Frequency).map_err(|m| s.err("frequency_unit", m))?;
            need_unit(s)?;
            let eps_inf = s.number("eps_inf")?.unwrap_or(1.0);
            let plasma = s
                .quantity("plasma_frequency", Dim::Frequency, Some(&fu))?
                .ok_or_else(|| s.err("plasma_frequency", "is missing; give the plasma frequency (the energy or frequency scale of the oscillator strengths)"))?;
            if plasma <= 0.0 {
                return Err(s.err("plasma_frequency", "must be positive"));
            }
            let drude = if name == "lorentz" {
                if s.has("damping") || s.has("drude_weight") {
                    return Err(s.err("damping", "belongs to a Drude term; use model = \"drude\" or \"drude_lorentz\""));
                }
                None
            } else {
                let g = s.quantity("damping", Dim::Frequency, Some(&fu))?.ok_or_else(|| s.err("damping", "is missing; a Drude term needs its damping (collision frequency)"))?;
                if g <= 0.0 {
                    return Err(s.err("damping", "must be positive (zero damping gives a lossless metal and a pole at zero frequency)"));
                }
                let f0 = s.number("drude_weight")?.unwrap_or(1.0);
                Some([f0, g])
            };
            let lorentz = if name == "drude" {
                if s.has("oscillators") {
                    return Err(s.err("oscillators", "belongs to Lorentz terms; use model = \"lorentz\" or \"drude_lorentz\""));
                }
                Vec::new()
            } else {
                let rows = match s.value("oscillators") {
                    Some(toml::Value::Array(r)) if !r.is_empty() => r,
                    _ => return Err(s.err("oscillators", "needs oscillators = [[f, resonance, damping], ...] with one row per Lorentz term (f is the weight, the other two are in frequency_unit)")),
                };
                let mut v = Vec::new();
                for (i, row) in rows.iter().enumerate() {
                    let toml::Value::Array(r) = row else { return Err(s.err("oscillators", format!("row {} must be [f, resonance, damping]", i + 1))) };
                    if r.len() != 3 {
                        return Err(s.err("oscillators", format!("row {} has {} entries; expected [f, resonance, damping]", i + 1, r.len())));
                    }
                    let x: Vec<f64> = r.iter().map(|e| s.number_of("oscillators", e)).collect::<Result<_, _>>()?;
                    if x[1] <= 0.0 || x[2] < 0.0 {
                        return Err(s.err("oscillators", format!("row {}: the resonance must be positive and the damping not negative", i + 1)));
                    }
                    v.push([x[0], x[1] * scale, x[2] * scale]);
                }
                v
            };
            Dispersion::Oscillators { eps_inf, plasma, drude, lorentz, frequency_unit: fu }
        }
        other => return Err(s.err("model", format!("unknown model '{other}'; use constant, sellmeier, cauchy, schott, herzberger, rii_formula_4, rii_formula_5, rii_formula_6, drude, lorentz, drude_lorentz or tabulated"))),
    };
    Ok(AxisModel {
        dispersion,
        source: Vec::new(),
        dn_dt,
        k_data,
        dn_dt_einstein_theta,
        dn_dt_range,
    })
}

fn fl(v: &[f64]) -> toml::Value {
    toml::Value::Array(v.iter().map(|x| toml::Value::Float(*x)).collect())
}

fn write_model(t: &mut toml::Table, m: &AxisModel) {
    t.insert("model".into(), m.dispersion.name().into());
    match &m.dispersion {
        Dispersion::Constant { n, k } => {
            t.insert("n".into(), toml::Value::Float(*n));
            if *k != 0.0 {
                t.insert("k".into(), toml::Value::Float(*k));
            }
        }
        Dispersion::Sellmeier { a, b, c } => {
            if *a != 1.0 {
                t.insert("A".into(), toml::Value::Float(*a));
            }
            t.insert("B".into(), fl(b));
            t.insert("C".into(), fl(c));
        }
        Dispersion::Cauchy { coefficients } => {
            t.insert("coefficients".into(), fl(coefficients));
        }
        Dispersion::Tabulated { lambda, n, k, interpolation, extrapolate } => {
            let with_k = k.iter().any(|x| *x != 0.0);
            let rows = lambda
                .iter()
                .enumerate()
                .map(|(i, l)| if with_k { fl(&[*l, n[i], k[i]]) } else { fl(&[*l, n[i]]) })
                .collect();
            t.insert("data".into(), toml::Value::Array(rows));
            if *interpolation == Interp::Linear {
                t.insert("interpolation".into(), "linear".into());
            }
            if *extrapolate == Extrapolate::Linear {
                t.insert("extrapolate".into(), "linear".into());
            }
        }
        Dispersion::Schott { coefficients }
        | Dispersion::Herzberger { coefficients }
        | Dispersion::RiiFormula4 { coefficients }
        | Dispersion::RiiFormula5 { coefficients }
        | Dispersion::RiiFormula6 { coefficients } => {
            t.insert("coefficients".into(), fl(coefficients));
        }
        Dispersion::Oscillators { eps_inf, plasma, drude, lorentz, frequency_unit } => {
            let s = freq_scale(frequency_unit);
            t.insert("frequency_unit".into(), frequency_unit.clone().into());
            if *eps_inf != 1.0 {
                t.insert("eps_inf".into(), toml::Value::Float(*eps_inf));
            }
            t.insert("plasma_frequency".into(), toml::Value::Float(plasma / s));
            if let Some([f0, g0]) = drude {
                t.insert("damping".into(), toml::Value::Float(g0 / s));
                if *f0 != 1.0 {
                    t.insert("drude_weight".into(), toml::Value::Float(*f0));
                }
            }
            if !lorentz.is_empty() {
                let rows = lorentz.iter().map(|[f, w0, g]| fl(&[*f, w0 / s, g / s])).collect();
                t.insert("oscillators".into(), toml::Value::Array(rows));
            }
        }
    }
    if let Some(v) = write_sources(&m.source) {
        t.insert("source".into(), v);
    }
    if m.dn_dt != 0.0 {
        t.insert("dn_dT".into(), format!("{} 1/K", m.dn_dt).into());
    }
    if let Some(theta) = m.dn_dt_einstein_theta {
        t.insert("dn_dT_einstein_theta".into(), format!("{theta} K").into());
    }
    if let Some([lo, hi]) = m.dn_dt_range {
        t.insert("dn_dT_range".into(), fl(&[lo, hi]));
    }
    if let Some((kl, kv)) = &m.k_data {
        let rows = kl.iter().enumerate().map(|(i, l)| fl(&[*l, kv[i]])).collect();
        t.insert("k_data".into(), toml::Value::Array(rows));
    }
}

/// Keys of the `[optical]` table itself (the model keys of an isotropic material live there too).
pub const OPTICAL_KEYS: &[&str] = &["wavelength_unit", "range", "reference_temperature", "ordinary", "extraordinary", "x", "y", "z"];

fn freq_scale(unit: &str) -> f64 {
    units::to_si(1.0, unit, Dim::Frequency).unwrap_or(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn silica() -> Optical {
        Optical {
            source: Vec::new(),
            wavelength_unit: "um".into(),
            unit_m: 1e-6,
            range: Some((0.21, 6.7)),
            reference_temperature: 293.15,
            axes: Axes::Isotropic(AxisModel {
                source: Vec::new(),
                dispersion: Dispersion::Sellmeier { a: 1.0, b: vec![0.6961663, 0.4079426, 0.8974794], c: vec![0.0684043_f64.powi(2), 0.1162414_f64.powi(2), 9.896161_f64.powi(2)] },
                dn_dt: 0.0,
                k_data: None,
                dn_dt_einstein_theta: None,
                dn_dt_range: None,
            }),
        }
    }

    #[test]
    fn malitson_silica_at_sodium_d_and_1550() {
        let o = silica();
        let n_d = o.nk(0, 0.5875618e-6, None).unwrap().0;
        assert!((n_d - 1.458464).abs() < 2e-6, "{n_d}");
        let n_1550 = o.nk(0, 1.55e-6, None).unwrap().0;
        assert!((n_1550 - 1.44402).abs() < 1e-4, "{n_1550}");
    }

    #[test]
    fn out_of_range_is_an_error_not_an_extrapolation() {
        let e = silica().nk(0, 0.1e-6, None).unwrap_err().to_string();
        assert!(e.contains("outside the range") && e.contains("0.21"), "{e}");
    }

    #[test]
    fn analytic_derivative_matches_finite_difference() {
        let o = silica();
        let (l, h) = (1.3e-6, 1e-12);
        let d = o.nk(0, l, None).unwrap().2;
        let fd = (o.nk(0, l + h, None).unwrap().0 - o.nk(0, l - h, None).unwrap().0) / (2.0 * h);
        assert!((d - fd).abs() < 1e-4 * fd.abs(), "{d} {fd}");
        // group index of silica at 1.55 um is about 1.4627
        let ng = o.group_index(0, 1.55e-6, None).unwrap();
        assert!((ng - 1.4627).abs() < 5e-4, "{ng}");
    }

    #[test]
    fn pchip_table_passes_through_points_and_stays_monotone() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [3.5, 3.4, 3.39, 3.2, 3.19];
        for (xi, yi) in x.iter().zip(y) {
            let (v, _) = interpolate(&x, &y, Interp::Pchip, *xi, false).unwrap();
            assert!((v - yi).abs() < 1e-12);
        }
        let mut prev = f64::INFINITY;
        for i in 0..=400 {
            let v = interpolate(&x, &y, Interp::Pchip, 1.0 + i as f64 * 0.01, false).unwrap().0;
            assert!(v <= prev + 1e-12);
            prev = v;
        }
        assert!(interpolate(&x, &y, Interp::Pchip, 5.1, false).is_none());
    }

    #[test]
    fn rii_formulas_4_5_6_match_finite_difference_derivatives() {
        // Formula 4 (SiO2/LiNbO3 form with negative IR pole power term):
        let f4 = Optical {
            source: Vec::new(),
            wavelength_unit: "um".into(),
            unit_m: 1e-6,
            range: Some((0.4, 5.0)),
            reference_temperature: 293.15,
            axes: Axes::Isotropic(AxisModel {
                dispersion: Dispersion::RiiFormula4 {
                    coefficients: vec![4.9048, 0.11768, 0.0, 0.0475, 1.0, 0.0, 0.0, 1.0, 1.0, -0.027169, 2.0],
                },
                source: Vec::new(),
                dn_dt: 0.0,
                k_data: None,
                dn_dt_einstein_theta: None,
                dn_dt_range: None,
            }),
        };
        let f5 = Optical {
            axes: Axes::Isotropic(AxisModel {
                dispersion: Dispersion::RiiFormula5 { coefficients: vec![2.391, 0.012, -2.0, 0.0005, -4.0] },
                source: Vec::new(),
                dn_dt: 0.0,
                k_data: None,
                dn_dt_einstein_theta: None,
                dn_dt_range: None,
            }),
            ..f4.clone()
        };
        let f6 = Optical {
            axes: Axes::Isotropic(AxisModel {
                dispersion: Dispersion::RiiFormula6 { coefficients: vec![0.0, 0.05792105, 238.0185, 0.00167917, 57.362] },
                source: Vec::new(),
                dn_dt: 0.0,
                k_data: None,
                dn_dt_einstein_theta: None,
                dn_dt_range: None,
            }),
            ..f4.clone()
        };
        for o in [f4, f5, f6] {
            let l = 1.55e-6;
            let h = 1e-11;
            let (_, _, dn, _) = o.nk(0, l, None).unwrap();
            let fd = (o.nk(0, l + h, None).unwrap().0 - o.nk(0, l - h, None).unwrap().0) / (2.0 * h);
            assert!((dn - fd).abs() < 1e-3 * fd.abs().max(1.0), "derivative mismatch: {dn} vs {fd}");
        }
    }

    #[test]
    fn k_data_and_einstein_dn_dt_work_with_sellmeier() {
        let mut o = silica();
        if let Axes::Isotropic(ref mut m) = o.axes {
            m.dn_dt = 1.0e-5;
            m.dn_dt_einstein_theta = Some(400.0);
            m.k_data = Some((vec![0.3, 0.5, 1.0, 1.5], vec![1e-3, 1e-4, 1e-6, 0.0]));
        }
        // Analytical n preserved alongside tabulated k
        let (n_500, k_500, _, _) = o.nk(0, 0.5e-6, None).unwrap();
        assert!((k_500 - 1e-4).abs() < 1e-12);
        assert!(n_500 > 1.45 && n_500 < 1.47);
        // Transparent tail past k_data returns k = 0 when within optical range
        let (_, k_1800, _, _) = o.nk(0, 1.8e-6, None).unwrap();
        assert_eq!(k_1800, 0.0);
        // Cryo Einstein saturation: slope -> 0 near 0 K
        let n_1k = o.nk(0, 1.55e-6, Some(1.0)).unwrap().0;
        let n_10k = o.nk(0, 1.55e-6, Some(10.0)).unwrap().0;
        let n_293 = o.nk(0, 1.55e-6, Some(293.15)).unwrap().0;
        assert!((n_10k - n_1k).abs() < 1e-18);
        assert!(n_10k < n_293);
    }
}
