// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Zero-thickness 2D material surface conductivity models (`[sheet_2d]` TOML v2 section, ADR 137):
//!
//! 1. **Analytical Kubo graphene (`Sheet2dModel::KuboGraphene`)**:
//!    - Intraband Drude-Boltzmann conductivity with exact finite-temperature chemical potential factor
//!      $\Delta_{\text{intra}}(\mu_c, T) = 2 k_B T \ln(2\cosh(\mu_c / (2 k_B T)))$ and optional
//!      out-of-plane magnetic field bias $B_0$ producing magneto-optic Dirac cyclotron Hall
//!      off-diagonals $\sigma_{xy} = -\sigma_{yx}$ ($\omega_c = e B_0 v_F^2 / \Delta_{\text{intra}}$).
//!    - Interband Kubo conductivity using the exact closed-form complex logarithm at $T \to 0$
//!      and nonsingular subtracted 32-point Gauss-Legendre quadrature with exact analytic tail
//!      integration at finite $T > 0$ (Falkovsky & Varlamov 2007; Gusynin, Sharapov & Carbotte 2007;
//!      Hanson 2008).
//! 2. **Anisotropic 2D sheets (`Sheet2dModel::Anisotropic`)**:
//!    - Principal-axis 2D surface conductivities $\sigma_{uu}(\omega), \sigma_{vv}(\omega)$ via
//!      constant baseline + 2D Drude + 2D Lorentz oscillator terms (e.g. black phosphorus, $\text{MoS}_2$,
//!      $\text{WSe}_2$, hexagonal boron nitride hyperbolic phonon-polariton sheets), rotated in-plane
//!      by $\theta_{2\text{D}}$:
//!      $$\boldsymbol{\sigma}_{2\text{D}}(\theta_{2\text{D}}) = \mathbf{R}(\theta_{2\text{D}})
//!      \operatorname{diag}(\sigma_{uu}, \sigma_{vv})\,\mathbf{R}(\theta_{2\text{D}})^T.$$
//! 3. **Tabulated 2D sheets (`Sheet2dModel::Tabulated`)**:
//!    - Shape-preserving Fritsch-Carlson PCHIP interpolation of tabulated complex surface
//!      conductivities $\sigma_{uu}(\omega), \sigma_{vv}(\omega)$ with literature provenance (`Source`)
//!      and valid temperature window (`T_range`).
//! 4. **Electrostatic gate-bias carrier-density mapping**:
//!    - Converts local electrostatic gate potential $V(\mathbf{r})$ (from `nyx-dfs` or uniform bias)
//!      and gate oxide capacitance per unit area $C_{\text{ox}} = \varepsilon_0 \varepsilon_{\text{ox}} / d_{\text{ox}}$
//!      into 2D carrier sheet density $n_s(\mathbf{r}) = C_{\text{ox}}(V(\mathbf{r}) - V_{\text{CNP}}) / e$
//!      and Dirac-cone chemical potential $\mu_c(\mathbf{r}) = \hbar v_F \sqrt{\pi |n_s(\mathbf{r})|}$.
//!
//! ### Time-Harmonic Sign Convention
//! Consistent with `nyx-ems` and `femcore::material::bcs`, all complex surface conductivities use the
//! $\exp(+i\omega t)$ phasor convention:
//! $$\boldsymbol{\sigma}_{2\text{D}}(\omega) = \boldsymbol{\sigma}'(\omega) - i\,\boldsymbol{\sigma}''(\omega),$$
//! where passive resistive dissipation satisfies $\text{Re}(\sigma_{uu}), \text{Re}(\sigma_{vv}) \ge 0$
//! and inductive (plasmonic/intraband) response has $\sigma'' > 0$ ($\text{Im}(\sigma) < 0$).

use super::bcs::{EV_TO_J, HBAR, K_B, MU0};
use super::reader::{MaterialError, Section};
use super::source::{read_sources, write_sources, Source};
use super::state_t::{clamp_temp_range, CheckedTempEval};
use super::units::{self, Dim};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Vacuum permittivity $\varepsilon_0$ in $\text{F/m}$.
pub const EPS0: f64 = 8.854_187_812_8e-12;
/// Speed of light in vacuum $c_0$ in $\text{m/s}$.
pub const C0: f64 = 299_792_458.0;
/// Free-space wave impedance $\eta_0 = \sqrt{\mu_0 / \varepsilon_0}$ in $\Omega$.
pub const ETA0: f64 = 376.730_313_668;
/// Elementary charge $e$ in Coulombs.
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Universal high-frequency interband optical sheet conductivity of monolayer graphene
/// $\sigma_0 = e^2 / (4\hbar) \approx 6.0853 \times 10^{-5}\text{ S}$.
pub const GRAPHENE_SIGMA0: f64 = E_CHARGE * E_CHARGE / (4.0 * HBAR);
/// Default monolayer graphene Dirac Fermi velocity $v_F = 1.0 \times 10^6\text{ m/s}$.
pub const DEFAULT_GRAPHENE_VF: f64 = 1.0e6;

/// Allowed TOML keys in the `[sheet_2d]` table.
pub const SHEET_2D_KEYS: &[&str] = &[
    "source",
    "model",
    "mu_c",
    "mu_c_unit",
    "tau",
    "tau_unit",
    "gamma",
    "gamma_unit",
    "scattering_rate",
    "scattering_rate_unit",
    "temperature",
    "temperature_unit",
    "T_range",
    "v_f",
    "v_f_unit",
    "b_0",
    "b_0_unit",
    "n_layers",
    "theta_deg",
    "theta_rad",
    "sigma_uu_re",
    "sigma_uu_im",
    "sigma_vv_re",
    "sigma_vv_im",
    "sigma_unit",
    "frequency_unit",
    "poles_u",
    "poles_v",
    "omega",
    "data",
    "data_uu",
    "data_vv",
];

// 32-point Gauss-Legendre nodes and weights on [-1, 1] (deterministic fixed-order quadrature).
const GL32_X: [f64; 32] = [
    -0.997_263_861_849_481_6,
    -0.985_611_511_545_268_4,
    -0.964_762_255_587_506_4,
    -0.934_906_075_937_739_7,
    -0.896_321_155_766_052_1,
    -0.849_367_613_732_57,
    -0.794_483_795_967_942_4,
    -0.732_182_118_740_289_7,
    -0.663_044_266_930_215_2,
    -0.587_715_757_240_762_3,
    -0.506_899_908_932_229_4,
    -0.421_351_276_130_635_33,
    -0.331_868_602_282_127_67,
    -0.239_287_362_252_137_06,
    -0.144_471_961_582_796_5,
    -0.048_307_665_687_738_32,
    0.048_307_665_687_738_32,
    0.144_471_961_582_796_5,
    0.239_287_362_252_137_06,
    0.331_868_602_282_127_67,
    0.421_351_276_130_635_33,
    0.506_899_908_932_229_4,
    0.587_715_757_240_762_3,
    0.663_044_266_930_215_2,
    0.732_182_118_740_289_7,
    0.794_483_795_967_942_4,
    0.849_367_613_732_57,
    0.896_321_155_766_052_1,
    0.934_906_075_937_739_7,
    0.964_762_255_587_506_4,
    0.985_611_511_545_268_4,
    0.997_263_861_849_481_6,
];

const GL32_W: [f64; 32] = [
    0.007_018_610_009_470_096,
    0.016_274_394_730_905_67,
    0.025_392_065_309_262_06,
    0.034_273_862_913_021_43,
    0.042_835_898_022_226_68,
    0.050_998_059_262_376_176,
    0.058_684_093_478_535_55,
    0.065_822_222_776_361_85,
    0.072_345_794_108_848_52,
    0.078_193_895_787_070_31,
    0.083_311_924_226_946_75,
    0.087_652_093_004_403_81,
    0.091_173_878_695_763_88,
    0.093_844_399_080_804_57,
    0.095_638_720_079_274_86,
    0.096_540_088_514_727_8,
    0.096_540_088_514_727_8,
    0.095_638_720_079_274_86,
    0.093_844_399_080_804_57,
    0.091_173_878_695_763_88,
    0.087_652_093_004_403_81,
    0.083_311_924_226_946_75,
    0.078_193_895_787_070_31,
    0.072_345_794_108_848_52,
    0.065_822_222_776_361_85,
    0.058_684_093_478_535_55,
    0.050_998_059_262_376_176,
    0.042_835_898_022_226_68,
    0.034_273_862_913_021_43,
    0.025_392_065_309_262_06,
    0.016_274_394_730_905_67,
    0.007_018_610_009_470_096,
];

fn gl32_integrate_c<F: Fn(f64) -> Complex64>(a: f64, b: f64, f: F) -> Complex64 {
    if b <= a {
        return Complex64::new(0.0, 0.0);
    }
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    let mut sum = Complex64::new(0.0, 0.0);
    for i in 0..32 {
        let x = mid + half * GL32_X[i];
        sum += f(x) * GL32_W[i];
    }
    sum * half
}

/// Shape-preserving Fritsch-Carlson PCHIP interpolation on strictly increasing `xs`.
/// Clamps to endpoint values outside `[xs[0], xs[n - 1]]`.
fn pchip_eval_sheet(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    debug_assert!(n >= 2 && ys.len() == n);
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    if n == 2 {
        let t = (x - xs[0]) / (xs[1] - xs[0]);
        return ys[0] + t * (ys[1] - ys[0]);
    }
    let mut h = Vec::with_capacity(n - 1);
    let mut delta = Vec::with_capacity(n - 1);
    for k in 0..n - 1 {
        let hk = xs[k + 1] - xs[k];
        h.push(hk);
        delta.push((ys[k + 1] - ys[k]) / hk);
    }
    let mut d = vec![0.0; n];
    let end_slope = |h0: f64, h1: f64, d0: f64, d1: f64| -> f64 {
        let s = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if s * d0 <= 0.0 {
            0.0
        } else if d0 * d1 <= 0.0 && s.abs() > 3.0 * d0.abs() {
            3.0 * d0
        } else {
            s
        }
    };
    d[0] = end_slope(h[0], h[1], delta[0], delta[1]);
    d[n - 1] = end_slope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
    for k in 1..n - 1 {
        if delta[k - 1] * delta[k] > 0.0 {
            let w1 = 2.0 * h[k] + h[k - 1];
            let w2 = h[k] + 2.0 * h[k - 1];
            d[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
        }
    }
    let mut lo = 0;
    let mut hi = n - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xs[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let hk = h[lo];
    let t = (x - xs[lo]) / hk;
    let t2 = t * t;
    let t3 = t2 * t;
    (2.0 * t3 - 3.0 * t2 + 1.0) * ys[lo]
        + (t3 - 2.0 * t2 + t) * hk * d[lo]
        + (-2.0 * t3 + 3.0 * t2) * ys[lo + 1]
        + (t3 - t2) * hk * d[lo + 1]
}

/// 2D Lorentz/Drude oscillator pole for principal-axis sheet conductivity:
/// $$\sigma_{\text{pole}}(\omega) = \frac{i \omega\, \sigma_{\text{weight}}}{\omega_0^2 - \omega^2 + i \gamma \omega}$$
/// (when $\omega_0 = 0$, reduces to a 2D Drude term $\frac{\sigma_{\text{weight}}}{\gamma + i\omega}$).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Sheet2dPole {
    /// Resonance angular frequency $\omega_0$ in $\text{rad/s}$ ($0$ for a pure 2D Drude pole).
    pub omega_0: f64,
    /// Damping rate $\gamma$ in $\text{rad/s}$ ($> 0$).
    pub gamma: f64,
    /// Oscillator spectral weight $S$ in $\text{S}\cdot\text{rad/s}$ (so at $\omega_0 = 0$, DC sheet conductivity is $S / \gamma$).
    pub weight: f64,
}

impl Sheet2dPole {
    /// Evaluates the oscillator contribution at angular frequency `omega` ($\text{rad/s}$) in $\exp(+i\omega t)$ convention.
    #[inline]
    pub fn eval(&self, omega: f64) -> Complex64 {
        let num = Complex64::new(0.0, omega * self.weight);
        let den = Complex64::new(self.omega_0 * self.omega_0 - omega * omega, self.gamma * omega);
        if den.norm_sqr() <= 1e-300 {
            Complex64::new(0.0, 0.0)
        } else {
            num / den
        }
    }
}

/// Principal-axis 2D surface conductivity model variants.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Sheet2dModel {
    /// Analytical Kubo monolayer (or N-layer) graphene model with intraband Drude + interband Kubo
    /// and optional out-of-plane magnetic field bias $B_0$ (Hall effect).
    KuboGraphene {
        /// Chemical potential (Fermi level) $\mu_c$ in Joules.
        mu_c: f64,
        /// Carrier momentum relaxation time $\tau$ in seconds ($\gamma = 1/\tau$).
        tau: f64,
        /// Reference temperature $T$ in Kelvin.
        temperature: f64,
        /// Dirac Fermi velocity $v_F$ in $\text{m/s}$ (default $1.0 \times 10^6\text{ m/s}$).
        v_f: f64,
        /// Out-of-plane static magnetic bias $B_0$ in Tesla (default $0.0\text{ T}$).
        b_0: f64,
        /// Effective number of decoupled graphene layers $N_{\text{layers}} \ge 1$ (default `1.0`).
        n_layers: f64,
    },
    /// Anisotropic 2D sheet with principal-axis constant baseline + Drude/Lorentz oscillator poles
    /// along local principal axes $(u, v)$ rotated in-plane by `theta_rad`.
    Anisotropic {
        /// Constant baseline complex sheet conductivity `[re, im]` $\sigma_{uu,0}$ in Siemens.
        sigma_uu_0: [f64; 2],
        /// Constant baseline complex sheet conductivity `[re, im]` $\sigma_{vv,0}$ in Siemens.
        sigma_vv_0: [f64; 2],
        /// Oscillator poles along principal axis $u$.
        poles_u: Vec<Sheet2dPole>,
        /// Oscillator poles along principal axis $v$.
        poles_v: Vec<Sheet2dPole>,
        /// Frequency unit string used in TOML representation (e.g. `"THz"` or `"rad/s"`).
        frequency_unit: String,
    },
    /// Tabulated complex 2D surface conductivity $\sigma_{uu}(\omega), \sigma_{vv}(\omega)$ interpolated via PCHIP.
    Tabulated {
        /// Strictly increasing angular frequency knots $\omega_k$ in $\text{rad/s}$.
        omega: Vec<f64>,
        /// Real part of $\sigma_{uu}(\omega_k)$ in Siemens ($\ge 0$).
        sigma_uu_re: Vec<f64>,
        /// Imaginary part of $\sigma_{uu}(\omega_k)$ in Siemens (in $\exp(+i\omega t)$ convention).
        sigma_uu_im: Vec<f64>,
        /// Real part of $\sigma_{vv}(\omega_k)$ in Siemens ($\ge 0$).
        sigma_vv_re: Vec<f64>,
        /// Imaginary part of $\sigma_{vv}(\omega_k)$ in Siemens (in $\exp(+i\omega t)$ convention).
        sigma_vv_im: Vec<f64>,
        /// Frequency unit string used in TOML representation.
        frequency_unit: String,
    },
}

/// Complete `[sheet_2d]` TOML v2 section describing a zero-thickness 2D material sheet.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Sheet2d {
    /// Underlying conductivity dispersion model.
    pub model: Sheet2dModel,
    /// In-plane rotation angle $\theta_{2\text{D}}$ of the principal axes $(u, v)$ in radians.
    pub theta_rad: f64,
    /// Optional valid temperature range `[T_min, T_max]` in Kelvin (ADR 133).
    pub t_range: Option<[f64; 2]>,
    /// Optional literature/measurement provenance records (ADR 132).
    pub source: Vec<Source>,
}

impl Sheet2d {
    /// Creates an analytical Kubo graphene sheet specification from SI parameters.
    pub fn kubo_graphene(
        mu_c_j: f64,
        tau_s: f64,
        temperature_k: f64,
        v_f: f64,
        b_0_t: f64,
    ) -> Result<Self, String> {
        if !mu_c_j.is_finite() {
            return Err("sheet_2d.mu_c must be finite".into());
        }
        if !tau_s.is_finite() || tau_s <= 0.0 {
            return Err("sheet_2d.tau must be finite and positive".into());
        }
        if !temperature_k.is_finite() || temperature_k < 0.0 {
            return Err("sheet_2d.temperature must be finite and non-negative".into());
        }
        if !v_f.is_finite() || v_f <= 0.0 {
            return Err("sheet_2d.v_f must be finite and positive".into());
        }
        if !b_0_t.is_finite() {
            return Err("sheet_2d.b_0 must be finite".into());
        }
        Ok(Self {
            model: Sheet2dModel::KuboGraphene {
                mu_c: mu_c_j,
                tau: tau_s,
                temperature: temperature_k,
                v_f,
                b_0: b_0_t,
                n_layers: 1.0,
            },
            theta_rad: 0.0,
            t_range: None,
            source: Vec::new(),
        })
    }

    /// Reads and validates a `[sheet_2d]` section from a TOML material file.
    pub fn read(s: &mut Section) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let model_str = s.string("model")?.unwrap_or_else(|| "kubo_graphene".to_string());

        let theta_deg = s.number("theta_deg")?;
        let theta_rad_raw = s.number("theta_rad")?;
        if theta_deg.is_some() && theta_rad_raw.is_some() {
            return Err(s.err("theta_rad", "specify either theta_deg or theta_rad, not both"));
        }
        let theta_rad = match (theta_deg, theta_rad_raw) {
            (Some(deg), None) => deg.to_radians(),
            (None, Some(rad)) => rad,
            _ => 0.0,
        };

        let t_range = match s.numbers("T_range")? {
            None => None,
            Some(v) if v.len() == 2 => {
                if v[0] < 0.0 || v[1] <= v[0] {
                    return Err(s.err(
                        "T_range",
                        format!("T_range must have 0 <= T_min < T_max (found [{}, {}])", v[0], v[1]),
                    ));
                }
                Some([v[0], v[1]])
            }
            Some(_) => return Err(s.err("T_range", "must be a two-element [T_min, T_max] array in K")),
        };

        let model = match model_str.as_str() {
            "kubo_graphene" | "kubo" => {
                let mu_c = s.quantity("mu_c", Dim::Energy, None)?.unwrap_or(0.3 * EV_TO_J);
                let tau_opt = s.quantity("tau", Dim::Time, None)?;
                let gamma_val = s.value("gamma").cloned();
                let _ = s.string("gamma_unit")?;
                let sr_val = s.value("scattering_rate").cloned();
                let _ = s.string("scattering_rate_unit")?;

                let tau = if let Some(t) = tau_opt {
                    if t <= 0.0 {
                        return Err(s.err("tau", format!("relaxation time tau must be positive (found {t})")));
                    }
                    t
                } else if gamma_val.is_some() || sr_val.is_some() {
                    let (key, val_opt, unit_key) = if let Some(v) = gamma_val {
                        ("gamma", v, "gamma_unit")
                    } else {
                        ("scattering_rate", sr_val.unwrap(), "scattering_rate_unit")
                    };
                    let decl_unit = s.string(unit_key)?;
                    let (num, unit) = match val_opt {
                        toml::Value::String(ref txt) => {
                            let (n, u) = units::split_quantity(txt).map_err(|m| s.err(key, m))?;
                            let eff_u = if u.is_empty() {
                                decl_unit.unwrap_or_default()
                            } else {
                                u
                            };
                            (n, eff_u)
                        }
                        other => {
                            let n = s.number_of(key, &other)?;
                            (n, decl_unit.unwrap_or_default())
                        }
                    };
                    if let Ok(e_j) = units::to_si(num, &unit, Dim::Energy) {
                        if e_j <= 0.0 {
                            return Err(s.err(key, "scattering rate must be positive"));
                        }
                        HBAR / (2.0 * e_j)
                    } else if let Ok(w_rads) = units::to_si(num, &unit, Dim::Frequency) {
                        if w_rads <= 0.0 {
                            return Err(s.err(key, "scattering rate must be positive"));
                        }
                        1.0 / w_rads
                    } else {
                        return Err(s.err(
                            key,
                            "must have energy units (e.g. \"1.3 meV\", Gamma = hbar/(2*tau)) or frequency units (e.g. \"2 THz\", gamma = 1/tau)",
                        ));
                    }
                } else {
                    0.5e-12
                };

                let temperature = s.quantity("temperature", Dim::Temperature, None)?.unwrap_or(300.0);
                if temperature < 0.0 {
                    return Err(s.err("temperature", "must be non-negative (K)"));
                }
                let v_f = s.quantity("v_f", Dim::Velocity, None)?.unwrap_or(DEFAULT_GRAPHENE_VF);
                if v_f <= 0.0 {
                    return Err(s.err("v_f", "Fermi velocity v_f must be positive"));
                }
                let b_0 = s.quantity("b_0", Dim::MagneticFluxDensity, None)?.unwrap_or(0.0);
                let n_layers = s.number("n_layers")?.unwrap_or(1.0);
                if n_layers <= 0.0 {
                    return Err(s.err("n_layers", "n_layers must be positive"));
                }
                Sheet2dModel::KuboGraphene {
                    mu_c,
                    tau,
                    temperature,
                    v_f,
                    b_0,
                    n_layers,
                }
            }
            "anisotropic" | "drude_lorentz_2d" | "constant" => {
                let sigma_unit = s.string("sigma_unit")?.unwrap_or_else(|| "S".to_string());
                let sigma_scale = units::to_si(1.0, &sigma_unit, Dim::SheetConductivity)
                    .map_err(|m| s.err("sigma_unit", m))?;
                let freq_unit = s.string("frequency_unit")?.unwrap_or_else(|| "rad/s".to_string());
                let freq_scale = units::to_si(1.0, &freq_unit, Dim::Frequency)
                    .map_err(|m| s.err("frequency_unit", m))?;

                let su_re = s
                    .quantity("sigma_uu_re", Dim::SheetConductivity, Some(&sigma_unit))?
                    .unwrap_or(0.0);
                if su_re < 0.0 {
                    return Err(s.err("sigma_uu_re", "real part of passive sheet conductivity must be >= 0"));
                }
                let su_im = s
                    .quantity("sigma_uu_im", Dim::SheetConductivity, Some(&sigma_unit))?
                    .unwrap_or(0.0);
                let sv_re = s
                    .quantity("sigma_vv_re", Dim::SheetConductivity, Some(&sigma_unit))?
                    .unwrap_or(su_re);
                if sv_re < 0.0 {
                    return Err(s.err("sigma_vv_re", "real part of passive sheet conductivity must be >= 0"));
                }
                let sv_im = s
                    .quantity("sigma_vv_im", Dim::SheetConductivity, Some(&sigma_unit))?
                    .unwrap_or(su_im);

                let read_poles = |sec: &mut Section, key: &str| -> Result<Vec<Sheet2dPole>, MaterialError> {
                    let Some(val) = sec.value(key) else {
                        return Ok(Vec::new());
                    };
                    let toml::Value::Array(rows) = val else {
                        return Err(sec.err(key, "must be a list of [omega_0, gamma, weight] rows"));
                    };
                    let mut poles = Vec::with_capacity(rows.len());
                    for (idx, row) in rows.iter().enumerate() {
                        let toml::Value::Array(cols) = row else {
                            return Err(sec.err(key, format!("row {idx} must be [omega_0, gamma, weight]")));
                        };
                        if cols.len() != 3 {
                            return Err(sec.err(key, format!("row {idx} must have 3 numbers [omega_0, gamma, weight]")));
                        }
                        let w0 = sec.number_of(key, &cols[0])? * freq_scale;
                        let gam = sec.number_of(key, &cols[1])? * freq_scale;
                        let wt = sec.number_of(key, &cols[2])? * sigma_scale * freq_scale;
                        if w0 < 0.0 {
                            return Err(sec.err(key, format!("row {idx}: omega_0 must be >= 0")));
                        }
                        if gam <= 0.0 {
                            return Err(sec.err(key, format!("row {idx}: gamma must be > 0")));
                        }
                        if wt < 0.0 {
                            return Err(sec.err(key, format!("row {idx}: weight must be >= 0 for passivity")));
                        }
                        poles.push(Sheet2dPole {
                            omega_0: w0,
                            gamma: gam,
                            weight: wt,
                        });
                    }
                    Ok(poles)
                };

                let poles_u = read_poles(s, "poles_u")?;
                let poles_v = if s.has("poles_v") {
                    read_poles(s, "poles_v")?
                } else if !s.has("sigma_vv_re") && !s.has("sigma_vv_im") {
                    poles_u.clone()
                } else {
                    Vec::new()
                };

                Sheet2dModel::Anisotropic {
                    sigma_uu_0: [su_re, su_im],
                    sigma_vv_0: [sv_re, sv_im],
                    poles_u,
                    poles_v,
                    frequency_unit: freq_unit,
                }
            }
            "tabulated" => {
                let sigma_unit = s.string("sigma_unit")?.unwrap_or_else(|| "S".to_string());
                let sigma_scale = units::to_si(1.0, &sigma_unit, Dim::SheetConductivity)
                    .map_err(|m| s.err("sigma_unit", m))?;
                let freq_unit = s.string("frequency_unit")?.unwrap_or_else(|| "THz".to_string());
                let freq_scale = units::to_si(1.0, &freq_unit, Dim::Frequency)
                    .map_err(|m| s.err("frequency_unit", m))?;

                let (omega, sigma_uu_re, sigma_uu_im, sigma_vv_re, sigma_vv_im) =
                    if let Some(val) = s.value("data") {
                        let toml::Value::Array(rows) = val else {
                            return Err(s.err("data", "must be an array of rows [freq, sigma_re, sigma_im] or [freq, su_re, su_im, sv_re, sv_im]"));
                        };
                        if rows.len() < 2 {
                            return Err(s.err("data", "tabulated sheet_2d requires at least 2 frequency points"));
                        }
                        let mut w = Vec::with_capacity(rows.len());
                        let mut u_re = Vec::with_capacity(rows.len());
                        let mut u_im = Vec::with_capacity(rows.len());
                        let mut v_re = Vec::with_capacity(rows.len());
                        let mut v_im = Vec::with_capacity(rows.len());
                        for (idx, r) in rows.iter().enumerate() {
                            let toml::Value::Array(cols) = r else {
                                return Err(s.err("data", format!("row {idx} must be a list of numbers")));
                            };
                            if cols.len() != 3 && cols.len() != 5 {
                                return Err(s.err("data", format!("row {idx} must have 3 (isotropic) or 5 (anisotropic) columns")));
                            }
                            let wk = s.number_of("data", &cols[0])? * freq_scale;
                            let ure = s.number_of("data", &cols[1])? * sigma_scale;
                            let uim = s.number_of("data", &cols[2])? * sigma_scale;
                            let (vre, vim) = if cols.len() == 5 {
                                (
                                    s.number_of("data", &cols[3])? * sigma_scale,
                                    s.number_of("data", &cols[4])? * sigma_scale,
                                )
                            } else {
                                (ure, uim)
                            };
                            w.push(wk);
                            u_re.push(ure);
                            u_im.push(uim);
                            v_re.push(vre);
                            v_im.push(vim);
                        }
                        (w, u_re, u_im, v_re, v_im)
                    } else if let Some(val_uu) = s.value("data_uu").cloned() {
                        type Col3 = (Vec<f64>, Vec<f64>, Vec<f64>);
                        let parse_3col = |sec: &Section, key: &str, v: &toml::Value| -> Result<Col3, MaterialError> {
                            let toml::Value::Array(rows) = v else {
                                return Err(sec.err(key, "must be an array of [freq, sigma_re, sigma_im] rows"));
                            };
                            if rows.len() < 2 {
                                return Err(sec.err(key, "requires at least 2 frequency points"));
                            }
                            let mut w = Vec::with_capacity(rows.len());
                            let mut re = Vec::with_capacity(rows.len());
                            let mut im = Vec::with_capacity(rows.len());
                            for (idx, r) in rows.iter().enumerate() {
                                let toml::Value::Array(cols) = r else {
                                    return Err(sec.err(key, format!("row {idx} must be [freq, sigma_re, sigma_im]")));
                                };
                                if cols.len() != 3 {
                                    return Err(sec.err(key, format!("row {idx} must have 3 columns [freq, sigma_re, sigma_im]")));
                                }
                                w.push(sec.number_of(key, &cols[0])? * freq_scale);
                                re.push(sec.number_of(key, &cols[1])? * sigma_scale);
                                im.push(sec.number_of(key, &cols[2])? * sigma_scale);
                            }
                            Ok((w, re, im))
                        };
                        let (w, u_re, u_im) = parse_3col(s, "data_uu", &val_uu)?;
                        let (v_re, v_im) = if let Some(val_vv) = s.value("data_vv").cloned() {
                            let (wv, vre, vim) = parse_3col(s, "data_vv", &val_vv)?;
                            if wv.len() != w.len()
                                || wv.iter().zip(&w).any(|(a, b)| (a - b).abs() > 1e-12 * b.abs().max(1.0))
                            {
                                return Err(s.err("data_vv", "data_vv frequency knots must match data_uu"));
                            }
                            (vre, vim)
                        } else {
                            (u_re.clone(), u_im.clone())
                        };
                        (w, u_re, u_im, v_re, v_im)
                    } else if let Some(w_raw) = s.numbers("omega")? {
                        if w_raw.len() < 2 {
                            return Err(s.err("omega", "tabulated sheet_2d requires at least 2 frequency points"));
                        }
                        let w: Vec<f64> = w_raw.into_iter().map(|f| f * freq_scale).collect();
                        let u_re: Vec<f64> = s
                            .numbers("sigma_uu_re")?
                            .ok_or_else(|| s.err("sigma_uu_re", "missing sigma_uu_re array for tabulated sheet_2d"))?
                            .into_iter()
                            .map(|v| v * sigma_scale)
                            .collect();
                        let u_im: Vec<f64> = s
                            .numbers("sigma_uu_im")?
                            .ok_or_else(|| s.err("sigma_uu_im", "missing sigma_uu_im array for tabulated sheet_2d"))?
                            .into_iter()
                            .map(|v| v * sigma_scale)
                            .collect();
                        let v_re: Vec<f64> = s
                            .numbers("sigma_vv_re")?
                            .map(|arr| arr.into_iter().map(|v| v * sigma_scale).collect())
                            .unwrap_or_else(|| u_re.clone());
                        let v_im: Vec<f64> = s
                            .numbers("sigma_vv_im")?
                            .map(|arr| arr.into_iter().map(|v| v * sigma_scale).collect())
                            .unwrap_or_else(|| u_im.clone());
                        if u_re.len() != w.len()
                            || u_im.len() != w.len()
                            || v_re.len() != w.len()
                            || v_im.len() != w.len()
                        {
                            return Err(s.err("omega", "all conductivity arrays must have the same length as omega"));
                        }
                        (w, u_re, u_im, v_re, v_im)
                    } else {
                        return Err(s.err("data", "tabulated sheet_2d requires `data` rows or `omega` + `sigma_uu_re` + `sigma_uu_im` arrays"));
                    };

                for k in 0..omega.len() {
                    if omega[k] <= 0.0 {
                        return Err(s.err("data", "tabulated frequencies must be strictly positive"));
                    }
                    if k > 0 && omega[k] <= omega[k - 1] {
                        return Err(s.err("data", "tabulated frequencies must be strictly increasing"));
                    }
                    if sigma_uu_re[k] < 0.0 || sigma_vv_re[k] < 0.0 {
                        return Err(s.err("data", "tabulated sheet conductivity real part must be >= 0 for passivity"));
                    }
                }

                Sheet2dModel::Tabulated {
                    omega,
                    sigma_uu_re,
                    sigma_uu_im,
                    sigma_vv_re,
                    sigma_vv_im,
                    frequency_unit: freq_unit,
                }
            }
            other => {
                return Err(s.err(
                    "model",
                    format!(
                        "unknown sheet_2d model '{other}'; expected \"kubo_graphene\", \"anisotropic\", or \"tabulated\""
                    ),
                ))
            }
        };

        s.finish(SHEET_2D_KEYS)?;
        Ok(Self {
            model,
            theta_rad,
            t_range,
            source,
        })
    }

    /// Serialises this `[sheet_2d]` specification into a TOML table (`format = 2`).
    pub fn write(&self, t: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            t.insert("source".into(), v);
        }
        if let Some([lo, hi]) = self.t_range {
            t.insert(
                "T_range".into(),
                toml::Value::Array(vec![toml::Value::Float(lo), toml::Value::Float(hi)]),
            );
        }
        if self.theta_rad != 0.0 {
            t.insert("theta_deg".into(), toml::Value::Float(self.theta_rad.to_degrees()));
        }
        match &self.model {
            Sheet2dModel::KuboGraphene {
                mu_c,
                tau,
                temperature,
                v_f,
                b_0,
                n_layers,
            } => {
                t.insert("model".into(), "kubo_graphene".into());
                t.insert("mu_c".into(), format!("{} eV", mu_c / EV_TO_J).into());
                t.insert("tau".into(), format!("{} ps", tau * 1e12).into());
                t.insert("temperature".into(), format!("{temperature} K").into());
                t.insert("v_f".into(), format!("{v_f} m/s").into());
                if *b_0 != 0.0 {
                    t.insert("b_0".into(), format!("{b_0} T").into());
                }
                if (*n_layers - 1.0).abs() > 1e-12 {
                    t.insert("n_layers".into(), toml::Value::Float(*n_layers));
                }
            }
            Sheet2dModel::Anisotropic {
                sigma_uu_0,
                sigma_vv_0,
                poles_u,
                poles_v,
                frequency_unit,
            } => {
                t.insert("model".into(), "anisotropic".into());
                t.insert("sigma_unit".into(), "S".into());
                t.insert("frequency_unit".into(), frequency_unit.clone().into());
                let f_scale = units::to_si(1.0, frequency_unit, Dim::Frequency).unwrap_or(1.0);
                t.insert("sigma_uu_re".into(), toml::Value::Float(sigma_uu_0[0]));
                t.insert("sigma_uu_im".into(), toml::Value::Float(sigma_uu_0[1]));
                t.insert("sigma_vv_re".into(), toml::Value::Float(sigma_vv_0[0]));
                t.insert("sigma_vv_im".into(), toml::Value::Float(sigma_vv_0[1]));
                if !poles_u.is_empty() {
                    let rows = poles_u
                        .iter()
                        .map(|p| {
                            toml::Value::Array(vec![
                                toml::Value::Float(p.omega_0 / f_scale),
                                toml::Value::Float(p.gamma / f_scale),
                                toml::Value::Float(p.weight / f_scale),
                            ])
                        })
                        .collect();
                    t.insert("poles_u".into(), toml::Value::Array(rows));
                }
                if !poles_v.is_empty() {
                    let rows = poles_v
                        .iter()
                        .map(|p| {
                            toml::Value::Array(vec![
                                toml::Value::Float(p.omega_0 / f_scale),
                                toml::Value::Float(p.gamma / f_scale),
                                toml::Value::Float(p.weight / f_scale),
                            ])
                        })
                        .collect();
                    t.insert("poles_v".into(), toml::Value::Array(rows));
                }
            }
            Sheet2dModel::Tabulated {
                omega,
                sigma_uu_re,
                sigma_uu_im,
                sigma_vv_re,
                sigma_vv_im,
                frequency_unit,
            } => {
                t.insert("model".into(), "tabulated".into());
                t.insert("sigma_unit".into(), "S".into());
                t.insert("frequency_unit".into(), frequency_unit.clone().into());
                let f_scale = units::to_si(1.0, frequency_unit, Dim::Frequency).unwrap_or(1.0);
                let rows = (0..omega.len())
                    .map(|k| {
                        toml::Value::Array(vec![
                            toml::Value::Float(omega[k] / f_scale),
                            toml::Value::Float(sigma_uu_re[k]),
                            toml::Value::Float(sigma_uu_im[k]),
                            toml::Value::Float(sigma_vv_re[k]),
                            toml::Value::Float(sigma_vv_im[k]),
                        ])
                    })
                    .collect();
                t.insert("data".into(), toml::Value::Array(rows));
            }
        }
    }

    /// Whether this 2D sheet exhibits a non-symmetric conductivity tensor ($\sigma_{xy} \neq \sigma_{yx}$),
    /// e.g. under a non-zero out-of-plane magnetic field bias $B_0 \neq 0$.
    pub fn is_non_symmetric(&self) -> bool {
        match &self.model {
            Sheet2dModel::KuboGraphene { b_0, .. } => b_0.abs() > 0.0,
            Sheet2dModel::Anisotropic { .. } | Sheet2dModel::Tabulated { .. } => false,
        }
    }

    /// Checks `temperature_k` against `t_range` (when present) and returns the effective temperature with clamping metadata.
    pub fn check_temperature(&self, temperature_k: f64) -> CheckedTempEval<f64> {
        match self.t_range {
            Some([lo, hi]) => {
                let (value, clamp) = clamp_temp_range(temperature_k, lo, hi);
                CheckedTempEval { value, clamp }
            }
            None => CheckedTempEval {
                value: temperature_k,
                clamp: None,
            },
        }
    }

    /// Evaluates the $2 \times 2$ tangential surface conductivity tensor
    /// $\boldsymbol{\sigma}_{2\text{D}}(\omega, T, \mu_c, B_0)$ in Siemens (`S`)
    /// at angular frequency `omega` ($\text{rad/s}$), optional temperature override `temp_k_override`,
    /// optional chemical potential override `mu_c_j_override` (e.g. from electrostatic gate bias),
    /// and optional magnetic field override `b_0_t_override`.
    pub fn eval_tensor(
        &self,
        omega: f64,
        temp_k_override: Option<f64>,
        mu_c_j_override: Option<f64>,
        b_0_t_override: Option<f64>,
    ) -> Result<[[Complex64; 2]; 2], String> {
        if !omega.is_finite() || omega <= 0.0 {
            return Err(format!(
                "sheet_2d evaluation angular frequency omega must be positive and finite (got {omega})"
            ));
        }
        if let Some(t_ov) = temp_k_override {
            if !t_ov.is_finite() || t_ov < 0.0 {
                return Err(format!(
                    "sheet_2d temperature override must be non-negative and finite (got {t_ov})"
                ));
            }
        }

        let (s_uu, s_vv, s_uv_hall) = match &self.model {
            Sheet2dModel::KuboGraphene {
                mu_c,
                tau,
                temperature,
                v_f,
                b_0,
                n_layers,
            } => {
                let mu_eff = mu_c_j_override.unwrap_or(*mu_c);
                let t_raw = temp_k_override.unwrap_or(*temperature);
                let t_eff = self.check_temperature(t_raw).value;
                let b_eff = b_0_t_override.unwrap_or(*b_0);
                let (s_xx, s_xy) = kubo_graphene_tensor(omega, t_eff, mu_eff, *tau, *v_f, b_eff)?;
                let scale = n_layers.max(1.0);
                (s_xx * scale, s_xx * scale, s_xy * scale)
            }
            Sheet2dModel::Anisotropic {
                sigma_uu_0,
                sigma_vv_0,
                poles_u,
                poles_v,
                ..
            } => {
                let mut su = Complex64::new(sigma_uu_0[0], sigma_uu_0[1]);
                for p in poles_u {
                    su += p.eval(omega);
                }
                let mut sv = Complex64::new(sigma_vv_0[0], sigma_vv_0[1]);
                for p in poles_v {
                    sv += p.eval(omega);
                }
                (su, sv, Complex64::new(0.0, 0.0))
            }
            Sheet2dModel::Tabulated {
                omega: w_grid,
                sigma_uu_re,
                sigma_uu_im,
                sigma_vv_re,
                sigma_vv_im,
                ..
            } => {
                let su = Complex64::new(
                    pchip_eval_sheet(w_grid, sigma_uu_re, omega),
                    pchip_eval_sheet(w_grid, sigma_uu_im, omega),
                );
                let sv = Complex64::new(
                    pchip_eval_sheet(w_grid, sigma_vv_re, omega),
                    pchip_eval_sheet(w_grid, sigma_vv_im, omega),
                );
                (su, sv, Complex64::new(0.0, 0.0))
            }
        };

        let (s, c) = self.theta_rad.sin_cos();
        let c2 = c * c;
        let s2 = s * s;
        let cs = c * s;
        // R(theta) [[s_uu, s_uv_hall], [-s_uv_hall, s_vv]] R(theta)^T
        // Note: an anti-symmetric Hall component [[0, h], [-h, 0]] is invariant under 2D SO(2) rotation!
        let s_xx = s_uu * c2 + s_vv * s2;
        let s_yy = s_uu * s2 + s_vv * c2;
        let s_sym_xy = (s_uu - s_vv) * cs;
        let t00 = s_xx;
        let t01 = s_sym_xy + s_uv_hall;
        let t10 = s_sym_xy - s_uv_hall;
        let t11 = s_yy;
        for row in [[t00, t01], [t10, t11]] {
            for val in row {
                if !val.re.is_finite() || !val.im.is_finite() {
                    return Err("non-finite 2D sheet conductivity tensor entry".into());
                }
            }
        }
        Ok([[t00, t01], [t10, t11]])
    }
}

/// Evaluates the analytical Kubo monolayer graphene conductivity components $(\sigma_{xx}, \sigma_{xy})$
/// in Siemens (`S`) at angular frequency `omega` ($\text{rad/s}$), temperature `temp_k` ($\text{K}$),
/// chemical potential `mu_c_j` ($\text{J}$), relaxation time `tau_s` ($\text{s}$), Fermi velocity `v_f` ($\text{m/s}$),
/// and out-of-plane magnetic field bias `b_0_t` ($\text{T}$), using the $\exp(+i\omega t)$ convention
/// ($\sigma = \sigma' - i\sigma''$, $\sigma' \ge 0$).
pub fn kubo_graphene_tensor(
    omega: f64,
    temp_k: f64,
    mu_c_j: f64,
    tau_s: f64,
    v_f: f64,
    b_0_t: f64,
) -> Result<(Complex64, Complex64), String> {
    if !omega.is_finite() || omega <= 0.0 {
        return Err(format!("omega must be finite and positive (got {omega})"));
    }
    if !temp_k.is_finite() || temp_k < 0.0 {
        return Err(format!("temperature must be finite and non-negative (got {temp_k})"));
    }
    if !mu_c_j.is_finite() {
        return Err(format!("mu_c must be finite (got {mu_c_j})"));
    }
    if !tau_s.is_finite() || tau_s <= 0.0 {
        return Err(format!("tau must be finite and positive (got {tau_s})"));
    }
    if !v_f.is_finite() || v_f <= 0.0 {
        return Err(format!("v_f must be finite and positive (got {v_f})"));
    }
    if !b_0_t.is_finite() {
        return Err(format!("b_0 must be finite (got {b_0_t})"));
    }

    let gamma = 1.0 / tau_s;
    let abs_mu = mu_c_j.abs();
    let kbt = K_B * temp_k;

    // Intraband effective chemical potential factor:
    // Delta_intra = 2 k_B T ln(2 cosh(mu_c / (2 k_B T))) = |mu_c| + 2 k_B T ln(1 + exp(-|mu_c| / (k_B T)))
    let delta_intra = if kbt <= 1e-30 * abs_mu.max(EV_TO_J * 1e-6) {
        abs_mu
    } else {
        let ratio = abs_mu / kbt;
        if ratio > 500.0 {
            abs_mu
        } else {
            abs_mu + 2.0 * kbt * (1.0 + (-ratio).exp()).ln()
        }
    };

    // Semiclassical / Dirac Kubo cyclotron frequency omega_c = e B_0 v_F^2 / Delta_intra
    let denom_energy = delta_intra.max(K_B * 1.0);
    let omega_c = E_CHARGE * b_0_t * v_f * v_f / denom_energy;

    let prefactor_intra = (E_CHARGE * E_CHARGE / (PI * HBAR * HBAR)) * delta_intra;
    // In exp(+i omega t), 1 / (gamma + i omega) gives Re > 0 and Im < 0 (inductive).
    let g_plus_iw = Complex64::new(gamma, omega);
    let denom_mag = g_plus_iw * g_plus_iw + Complex64::new(omega_c * omega_c, 0.0);
    let sigma_intra_xx = prefactor_intra * g_plus_iw / denom_mag;
    let sigma_intra_xy = -prefactor_intra * Complex64::new(omega_c, 0.0) / denom_mag;

    let sigma_inter = kubo_graphene_interband(omega, temp_k, abs_mu, gamma);
    let sigma_xx = sigma_intra_xx + sigma_inter;
    let sigma_xy = sigma_intra_xy;

    if !sigma_xx.re.is_finite()
        || !sigma_xx.im.is_finite()
        || !sigma_xy.re.is_finite()
        || !sigma_xy.im.is_finite()
    {
        return Err("non-finite Kubo graphene conductivity".into());
    }
    Ok((sigma_xx, sigma_xy))
}

/// Evaluates the isotropic ($B_0 = 0$) analytical Kubo monolayer graphene surface conductivity
/// $\sigma(\omega, T, \mu_c, \tau) = \sigma_{\text{intra}} + \sigma_{\text{inter}}$ in Siemens (`S`).
pub fn kubo_graphene_isotropic(
    omega: f64,
    temp_k: f64,
    mu_c_j: f64,
    tau_s: f64,
) -> Result<Complex64, String> {
    let (s_xx, _) = kubo_graphene_tensor(omega, temp_k, mu_c_j, tau_s, DEFAULT_GRAPHENE_VF, 0.0)?;
    Ok(s_xx)
}

/// Closed-form $T = 0$ and subtracted finite-$T$ Gauss-Legendre interband Kubo conductivity
/// in the $\exp(+i\omega t)$ convention ($\text{Re}(\sigma_{\text{inter}}) \ge 0$).
pub fn kubo_graphene_interband(omega: f64, temp_k: f64, abs_mu_j: f64, gamma: f64) -> Complex64 {
    let hw = HBAR * omega;
    let hgamma = HBAR * gamma;
    // w = hbar * (omega - i * gamma) in exp(+i omega t)
    let w = Complex64::new(hw, -hgamma);

    // Closed-form T = 0 expression:
    // sigma_inter,0 = -i * e^2 / (4 pi hbar) * ln((2|mu_c| - w) / (2|mu_c| + w))
    let kbt = K_B * temp_k;
    let scale_e = abs_mu_j.max(hw).max(EV_TO_J * 1e-6);
    if kbt <= 1e-4 * scale_e {
        let z_minus = Complex64::new(2.0 * abs_mu_j, 0.0) - w;
        let z_plus = Complex64::new(2.0 * abs_mu_j, 0.0) + w;
        let log_ratio = (z_minus / z_plus).ln();
        return Complex64::new(0.0, -GRAPHENE_SIGMA0 / PI) * log_ratio;
    }

    // Finite-temperature Fermi-Dirac interband population difference:
    // H(eps) = sinh(eps / kBT) / (cosh(mu_c / kBT) + cosh(eps / kBT))
    // Stable evaluation without overflow for arbitrarily large eps / kBT or mu_c / kBT:
    let h_fn = |eps: f64| -> f64 {
        let x = eps / kbt;
        let m = abs_mu_j / kbt;
        if x <= 0.0 {
            return 0.0;
        }
        let f_neg = 1.0 / (1.0 + (-(x + m)).exp());
        let f_pos = if x - m > 500.0 {
            0.0
        } else if x - m < -500.0 {
            1.0
        } else {
            1.0 / (1.0 + (x - m).exp())
        };
        (f_neg - f_pos).clamp(0.0, 1.0)
    };

    let eps_pole = 0.5 * hw;
    let h_pole = h_fn(eps_pole);
    let w2 = w * w;

    let integrand = |eps: f64| -> Complex64 {
        let dh = h_fn(eps) - h_pole;
        if dh.abs() <= 1e-16 && (eps - eps_pole).abs() <= 1e-12 * scale_e {
            return Complex64::new(0.0, 0.0);
        }
        let den = w2 - Complex64::new(4.0 * eps * eps, 0.0);
        Complex64::new(dh, 0.0) / den
    };

    // Upper cutoff where H(eps) == 1 to machine precision (eps >= |mu_c| + 45 k_B T and eps >= 2.5 * eps_pole)
    let eps_max = (abs_mu_j + 45.0 * kbt).max(2.5 * eps_pole + 20.0 * kbt);

    // Split [0, eps_max] at key transition landmarks
    let mut knots = vec![
        0.0,
        (abs_mu_j - 12.0 * kbt).max(0.0),
        (abs_mu_j - 3.0 * kbt).max(0.0),
        abs_mu_j,
        (abs_mu_j + 3.0 * kbt).min(eps_max),
        (abs_mu_j + 12.0 * kbt).min(eps_max),
        (eps_pole - 4.0 * hgamma).max(0.0).min(eps_max),
        eps_pole.min(eps_max),
        (eps_pole + 4.0 * hgamma).min(eps_max),
        eps_max,
    ];
    knots.sort_by(|a, b| a.total_cmp(b));

    let mut integral = Complex64::new(0.0, 0.0);
    for pair in knots.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b > a + 1e-15 * eps_max {
            integral += gl32_integrate_c(a, b, integrand);
        }
    }

    // Exact analytic tail on [eps_max, infinity) where H(eps) = 1:
    // \int_{eps_max}^\infty (1 - h_pole) / (w^2 - 4 eps^2) d(eps)
    //   = ((1 - h_pole) / (4 w)) * ln((2 eps_max - w) / (2 eps_max + w))
    let z_tail_minus = Complex64::new(2.0 * eps_max, 0.0) - w;
    let z_tail_plus = Complex64::new(2.0 * eps_max, 0.0) + w;
    let tail = Complex64::new(1.0 - h_pole, 0.0) * (z_tail_minus / z_tail_plus).ln() / (4.0 * w);
    integral += tail;

    // sigma_inter = (e^2 / 4 hbar) * [ H(hbar omega / 2) - (4 i w / pi) * integral ]
    let bracket = Complex64::new(h_pole, 0.0) - Complex64::new(0.0, 4.0 / PI) * w * integral;
    Complex64::new(GRAPHENE_SIGMA0, 0.0) * bracket
}

/// Computes gate oxide capacitance per unit area $C_{\text{ox}} = \varepsilon_0 \varepsilon_{\text{ox}} / d_{\text{ox}}$
/// in $\text{F/m}^2$ from relative oxide permittivity $\varepsilon_{\text{ox}} > 0$ and thickness $d_{\text{ox}} > 0$ ($\text{m}$).
pub fn gate_oxide_capacitance(eps_r_ox: f64, d_ox_m: f64) -> Result<f64, String> {
    if !eps_r_ox.is_finite() || eps_r_ox <= 0.0 {
        return Err(format!(
            "oxide relative permittivity eps_r_ox must be finite and positive (got {eps_r_ox})"
        ));
    }
    if !d_ox_m.is_finite() || d_ox_m <= 0.0 {
        return Err(format!(
            "oxide thickness d_ox_m must be finite and positive (got {d_ox_m})"
        ));
    }
    Ok(EPS0 * eps_r_ox / d_ox_m)
}

/// Maps electrostatic gate potential $V_{\text{gate}}$ ($\text{V}$), charge-neutrality offset $V_{\text{CNP}}$ ($\text{V}$),
/// gate capacitance per unit area $C_{\text{ox}}$ ($\text{F/m}^2$), and Dirac Fermi velocity $v_F$ ($\text{m/s}$)
/// to `(n_s_m2, mu_c_j)`:
/// $$n_s = \frac{C_{\text{ox}} (V_{\text{gate}} - V_{\text{CNP}})}{e}, \qquad
///   \mu_c = \hbar v_F \sqrt{\pi |n_s|}.$$
pub fn chemical_potential_from_gate_potential(
    v_gate: f64,
    v_cnp: f64,
    c_ox_f_per_m2: f64,
    v_f: f64,
) -> Result<(f64, f64), String> {
    if !v_gate.is_finite() || !v_cnp.is_finite() {
        return Err("gate potential and V_CNP must be finite".into());
    }
    if !c_ox_f_per_m2.is_finite() || c_ox_f_per_m2 < 0.0 {
        return Err(format!(
            "gate oxide capacitance C_ox must be finite and non-negative (got {c_ox_f_per_m2})"
        ));
    }
    if !v_f.is_finite() || v_f <= 0.0 {
        return Err(format!(
            "Fermi velocity v_f must be finite and positive (got {v_f})"
        ));
    }
    let n_s = c_ox_f_per_m2 * (v_gate - v_cnp) / E_CHARGE;
    let mu_c_j = HBAR * v_f * (PI * n_s.abs()).sqrt();
    Ok((n_s, mu_c_j))
}

/// Closed-form normal-incidence Fresnel transmission, reflection, and absorption $(T, R, A)$
/// of a zero-thickness isotropic conductive sheet with surface conductivity $\sigma(\omega)$ (in Siemens)
/// suspended at the interface between half-spaces of refractive indices $n_1, n_2 > 0$:
/// $$t = \frac{2 n_1}{n_1 + n_2 + \eta_0 \sigma}, \quad r = \frac{n_1 - n_2 - \eta_0 \sigma}{n_1 + n_2 + \eta_0 \sigma},$$
/// $$T = \frac{n_2}{n_1} |t|^2, \quad R = |r|^2, \quad A = 1 - T - R = \frac{4 n_1 \eta_0 \operatorname{Re}(\sigma)}{|n_1 + n_2 + \eta_0 \sigma|^2}.$$
pub fn sheet2d_fresnel_normal_tra(sigma: Complex64, n1: f64, n2: f64) -> (f64, f64, f64) {
    let eta0 = (MU0 / EPS0).sqrt();
    let denom = Complex64::new(n1 + n2, 0.0) + eta0 * sigma;
    let t_amp = Complex64::new(2.0 * n1, 0.0) / denom;
    let r_amp = (Complex64::new(n1 - n2, 0.0) - eta0 * sigma) / denom;
    let t_pow = (n2 / n1) * t_amp.norm_sqr();
    let r_pow = r_amp.norm_sqr();
    let a_pow = 1.0 - t_pow - r_pow;
    (t_pow, r_pow, a_pow)
}

/// Non-retarded (quasi-static, $q_{\text{spp}} \gg k_0$) 2D TM graphene surface plasmon dispersion
/// wavenumber $q_{\text{spp}}$ ($\text{rad/m}$) at a sheet with conductivity $\sigma(\omega)$
/// between dielectrics $\varepsilon_{r,1}, \varepsilon_{r,2}$ ($\bar{\varepsilon}_r = (\varepsilon_{r,1} + \varepsilon_{r,2})/2$)
/// in the $\exp(+i\omega t - i q_{\text{spp}} y)$ convention ($\operatorname{Re}(q_{\text{spp}}) > 0$, $\operatorname{Im}(q_{\text{spp}}) \le 0$):
/// $$q_{\text{spp}} = \frac{2\omega\varepsilon_0 \bar{\varepsilon}_r}{i\,\sigma(\omega)} = \frac{-i\,2\omega\varepsilon_0 \bar{\varepsilon}_r}{\sigma(\omega)}.$$
pub fn tm_plasmon_wavenumber_nonretarded(
    sigma: Complex64,
    omega: f64,
    eps_r1: f64,
    eps_r2: f64,
) -> Complex64 {
    let bar_eps_r = 0.5 * (eps_r1 + eps_r2);
    let num = Complex64::new(0.0, -2.0 * omega * EPS0 * bar_eps_r);
    num / sigma
}

/// Exact retarded 2D TM surface plasmon propagation constant $\beta_{\text{spp}}$ ($\text{rad/m}$)
/// for a conductive sheet embedded in a symmetric dielectric medium $\varepsilon_r$:
/// $$\beta_{\text{spp}} = k_0 \sqrt{\varepsilon_r - \left(\frac{2\varepsilon_r}{\eta_0 \sigma(\omega)}\right)^2},$$
/// chosen in the right half-plane $\operatorname{Re}(\beta_{\text{spp}}) \ge 0$.
pub fn tm_plasmon_wavenumber_retarded_symmetric(
    sigma: Complex64,
    omega: f64,
    eps_r: f64,
) -> Complex64 {
    let k0 = omega / C0;
    let eta0 = (MU0 / EPS0).sqrt();
    let ratio = Complex64::new(2.0 * eps_r, 0.0) / (eta0 * sigma);
    let mut beta = k0 * (Complex64::new(eps_r, 0.0) - ratio * ratio).sqrt();
    if beta.re < 0.0 {
        beta = -beta;
    }
    beta
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::MaterialSpec;

    #[test]
    fn kubo_graphene_intraband_and_universal_interband_limits() {
        // 1. At low frequency (1 THz) with mu_c = 0.4 eV, tau = 0.5 ps, T = 300 K:
        // intraband Drude dominates, Re(sigma) > 0, Im(sigma) < 0 (inductive plasmonic sheet).
        let mu_c = 0.4 * EV_TO_J;
        let tau = 0.5e-12;
        let omega_1thz = 2.0 * PI * 1.0e12;
        let s_1thz = kubo_graphene_isotropic(omega_1thz, 300.0, mu_c, tau).unwrap();
        assert!(s_1thz.re > 1.0e-3, "Re(sigma) at 1 THz = {}", s_1thz.re);
        assert!(s_1thz.im < -1.0e-3, "Im(sigma) at 1 THz = {}", s_1thz.im);

        // 2. At optical frequency hbar*omega = 2.5 eV >> 2*mu_c = 0.8 eV:
        // interband conductivity approaches the universal monolayer limit sigma_0 = e^2 / (4 hbar).
        let omega_opt = 2.5 * EV_TO_J / HBAR;
        let s_opt_0k = kubo_graphene_isotropic(omega_opt, 0.0, mu_c, tau).unwrap();
        let s_opt_300k = kubo_graphene_isotropic(omega_opt, 300.0, mu_c, tau).unwrap();
        assert!(
            (s_opt_0k.re - GRAPHENE_SIGMA0).abs() < 0.02 * GRAPHENE_SIGMA0,
            "0 K optical Re(sigma) = {}, want ~{}",
            s_opt_0k.re,
            GRAPHENE_SIGMA0
        );
        assert!(
            (s_opt_300k.re - GRAPHENE_SIGMA0).abs() < 0.02 * GRAPHENE_SIGMA0,
            "300 K optical Re(sigma) = {}, want ~{}",
            s_opt_300k.re,
            GRAPHENE_SIGMA0
        );

        // 3. Continuity between T = 0 closed-form and low-T (T = 10 K) subtracted quadrature:
        let omega_mid = 2.0 * PI * 30.0e12;
        let s_0 = kubo_graphene_isotropic(omega_mid, 0.0, mu_c, tau).unwrap();
        let s_10k = kubo_graphene_isotropic(omega_mid, 10.0, mu_c, tau).unwrap();
        assert!(
            (s_0 - s_10k).norm() < 5e-3 * s_0.norm(),
            "T=0 vs T=10K continuity: {s_0:?} vs {s_10k:?}"
        );
    }

    #[test]
    fn kubo_graphene_hall_bias_and_gate_potential_mapping() {
        let mu_c = 0.25 * EV_TO_J;
        let tau = 0.4e-12;
        let omega = 2.0 * PI * 5.0e12;
        let sheet = Sheet2d::kubo_graphene(mu_c, tau, 300.0, DEFAULT_GRAPHENE_VF, 5.0).unwrap();
        assert!(sheet.is_non_symmetric());
        let tensor = sheet.eval_tensor(omega, None, None, None).unwrap();
        assert_eq!(tensor[0][0], tensor[1][1]);
        assert!((tensor[0][1] + tensor[1][0]).norm() < 1e-18, "Hall anti-symmetry sigma_xy = -sigma_yx");
        assert!(tensor[0][1].norm() > 1e-5, "non-zero Hall off-diagonal at B_0 = 5 T");

        // Gate oxide mapping: 10 nm Al2O3 (eps_r = 9.0) at V_g - V_CNP = 2.0 V
        let c_ox = gate_oxide_capacitance(9.0, 10e-9).unwrap();
        let (n_s, mu_c_gated) =
            chemical_potential_from_gate_potential(2.0, 0.0, c_ox, DEFAULT_GRAPHENE_VF).unwrap();
        assert!(n_s > 0.0);
        let mu_ev = mu_c_gated / EV_TO_J;
        assert!(mu_ev > 0.35 && mu_ev < 0.45, "mu_c(2V) = {mu_ev} eV");
    }

    #[test]
    fn toml_parsing_and_roundtrip_for_all_sheet2d_models() {
        let toml_kubo = r#"
format = 2
name = "graphene_test"
kind = "semiconductor"
[sheet_2d]
model = "kubo_graphene"
mu_c = "0.35 eV"
tau = "0.5 ps"
temperature = "300 K"
v_f = "1.0e6 m/s"
b_0 = "1.5 T"
T_range = [4.0, 1000.0]
"#;
        let spec = MaterialSpec::from_toml(toml_kubo).unwrap();
        let sh = spec.sheet_2d.as_ref().unwrap();
        assert_eq!(sh.t_range, Some([4.0, 1000.0]));
        assert!(sh.is_non_symmetric());
        let round = MaterialSpec::from_toml(&spec.to_toml()).unwrap();
        assert_eq!(round.sheet_2d, spec.sheet_2d);

        // Anisotropic 2D sheet rotated by 30 deg
        let toml_aniso = r#"
format = 2
name = "bp_2d"
[sheet_2d]
model = "anisotropic"
theta_deg = 30.0
sigma_unit = "mS"
frequency_unit = "THz"
sigma_uu_re = 0.2
sigma_vv_re = 0.05
poles_u = [[0.0, 2.0, 4.0]]
poles_v = [[0.0, 2.0, 1.0]]
"#;
        let spec_a = MaterialSpec::from_toml(toml_aniso).unwrap();
        let sh_a = spec_a.sheet_2d.as_ref().unwrap();
        assert!(!sh_a.is_non_symmetric());
        let t_a = sh_a.eval_tensor(2.0 * PI * 10.0e12, None, None, None).unwrap();
        assert!((t_a[0][1] - t_a[1][0]).norm() < 1e-18, "rotated anisotropic sheet is symmetric");
        assert!(t_a[0][1].norm() > 1e-6, "30-deg rotation produces symmetric off-diagonal coupling");
    }
}

