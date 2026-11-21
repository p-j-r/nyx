// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Arbitrary-purity BCS complex conductivity and thin-film surface impedance evaluator
//! (Zimmermann et al., Physica C 183, 99-104, 1991; Mattis & Bardeen, Phys. Rev. 111, 412, 1958;
//! Park et al., npj Quantum Inf. 2026) extended with sub-gap mK loss mechanisms (`dynes_gamma`,
//! `r_res`, `x_qp_floor`, `nominal_thickness`, `film_model`).
//!
//! Time-harmonic convention: `exp(+i omega t)`, so the superconducting complex conductivity is
//! `sigma(omega, T) = sigma_1(omega, T) - i sigma_2(omega, T)` with `sigma_1 >= 0` and `sigma_2 > 0`
//! for `T < T_c`, and the sheet surface impedance is `Z_s = R_s + i omega L_k_sq` with `R_s >= 0`
//! and `L_k_sq > 0`.

use super::reader::{MaterialError, Section};
use super::source::{read_sources, write_sources, Source};
use super::units::Dim;
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Boltzmann constant `k_B` in J/K (SI exact).
pub const K_B: f64 = 1.380_649e-23;
/// Reduced Planck constant `hbar` in J*s.
pub const HBAR: f64 = 1.054_571_817e-34;
/// Elementary charge / eV-to-Joule conversion factor in J/eV (SI exact).
pub const EV_TO_J: f64 = 1.602_176_634e-19;
/// Vacuum magnetic permeability `mu_0` in H/m.
pub const MU0: f64 = 4.0 * PI * 1e-7;
/// Default weak-coupling BCS zero-temperature gap ratio `Delta(0) / (k_B T_c)`.
pub const BCS_GAP_RATIO: f64 = 1.7638;
/// Exact weak-coupling BCS universal constant `pi * exp(-gamma_E)` where `Delta(T -> T_c) -> 0`.
const BCS_EXACT_PI_EXP_GAMMA: f64 = 1.763_877_368_839_073_8;

/// Positive abscissas of the 32-point Gauss-Legendre quadrature rule on `[-1, 1]`.
#[allow(clippy::excessive_precision)]
const GL32_X: [f64; 16] = [
    0.048_307_665_687_738_316,
    0.144_471_961_582_796_493,
    0.239_287_362_252_137_075,
    0.331_868_602_282_127_650,
    0.421_351_276_130_635_345,
    0.506_899_908_932_229_390,
    0.587_715_757_240_762_329,
    0.663_044_266_930_215_201,
    0.732_182_118_740_289_680,
    0.794_483_795_967_942_407,
    0.849_367_613_732_569_970,
    0.896_321_155_766_052_124,
    0.934_906_075_937_739_689,
    0.964_762_255_587_506_431,
    0.985_611_511_545_268_335,
    0.997_263_861_849_481_564,
];

/// Weights corresponding to `GL32_X` for the 32-point Gauss-Legendre rule on `[-1, 1]`.
#[allow(clippy::excessive_precision)]
const GL32_W: [f64; 16] = [
    0.096_540_088_514_727_801,
    0.095_638_720_079_274_859,
    0.093_844_399_080_804_566,
    0.091_173_878_695_763_885,
    0.087_652_093_004_403_811,
    0.083_311_924_226_946_755,
    0.078_193_895_787_070_306,
    0.072_345_794_108_848_506,
    0.065_822_222_776_361_847,
    0.058_684_093_478_535_547,
    0.050_998_059_262_376_176,
    0.042_835_898_022_226_681,
    0.034_273_862_913_021_433,
    0.025_392_065_309_262_059,
    0.016_274_394_730_905_671,
    0.007_018_610_009_470_097,
];

/// Thin-film boundary model mapping bulk complex conductivity `sigma` and film thickness `d`
/// to effective sheet surface impedance `Z_s = R_s + i omega L_k_sq`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum FilmModel {
    /// Symmetric two-sided current distribution: `Z_s = (gamma / (2 sigma)) coth(gamma d / 2)`.
    TwoSided,
    /// Kerr single-sided / backed film formula (Park et al. 2026 Eq. 4): `Z_s = (gamma / sigma) coth(gamma d)`.
    #[default]
    Kerr,
    /// Thin-film local limit (`d << lambda_eff`): `Z_s = 1 / (sigma d)`.
    ThinLimit,
}

impl FilmModel {
    pub fn parse(s: &str) -> Option<FilmModel> {
        match s.trim().to_ascii_lowercase().as_str() {
            "two_sided" | "twosided" | "two-sided" => Some(FilmModel::TwoSided),
            "kerr" | "one_sided" | "onesided" => Some(FilmModel::Kerr),
            "thin_limit" | "thinlimit" | "thin-limit" | "thin" => Some(FilmModel::ThinLimit),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            FilmModel::TwoSided => "two_sided",
            FilmModel::Kerr => "kerr",
            FilmModel::ThinLimit => "thin_limit",
        }
    }
}

/// Superconducting normal-state, gap, mean-free-path, and sub-gap mK loss parameters (`[superconducting]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Superconducting {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<Source>,
    /// Critical temperature `T_c` in K.
    pub tc: f64,
    /// Zero-temperature superconducting gap `Delta(0)` in J (defaults to `1.7638 k_B T_c`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap_0: Option<f64>,
    /// Clean-limit London penetration depth `lambda_L` in m.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lambda_l: Option<f64>,
    /// Pippard/BCS coherence length `xi_0` in m.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub xi_0: Option<f64>,
    /// Electron mean free path `ell` in m.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mean_free_path: Option<f64>,
    /// Normal-state electrical conductivity `sigma_n` just above `T_c` in S/m.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sigma_n: Option<f64>,
    /// Dynes pair-breaking broadening `Gamma` in J (produces finite sub-gap density of states at mK).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dynes_gamma: Option<f64>,
    /// Non-BCS residual sheet surface resistance floor `R_res` in Ohm/sq.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r_res: Option<f64>,
    /// Non-equilibrium normalized quasiparticle density floor `x_qp = n_qp / (2 N_0 Delta)`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub x_qp_floor: Option<f64>,
    /// Default nominal film thickness in m.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nominal_thickness: Option<f64>,
    /// Default thin-film geometry model when evaluated at default thickness.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub film_model: Option<FilmModel>,
}

pub const SC_KEYS: &[&str] = &[
    "source",
    "tc",
    "tc_unit",
    "gap_0",
    "gap_0_unit",
    "lambda_L",
    "lambda_L_unit",
    "xi_0",
    "xi_0_unit",
    "mean_free_path",
    "mean_free_path_unit",
    "sigma_n",
    "sigma_n_unit",
    "dynes_gamma",
    "dynes_gamma_unit",
    "r_res",
    "r_res_unit",
    "x_qp_floor",
    "nominal_thickness",
    "nominal_thickness_unit",
    "film_model",
];

fn positive_opt(s: &Section, key: &str, v: Option<f64>, what: &str) -> Result<Option<f64>, MaterialError> {
    match v {
        Some(x) if !x.is_finite() || x <= 0.0 => Err(s.err(key, format!("{what} must be positive (found {x})"))),
        other => Ok(other),
    }
}

fn nonneg_opt(s: &Section, key: &str, v: Option<f64>, what: &str) -> Result<Option<f64>, MaterialError> {
    match v {
        Some(x) if !x.is_finite() || x < 0.0 => Err(s.err(key, format!("{what} must be non-negative (found {x})"))),
        other => Ok(other),
    }
}

impl Superconducting {
    /// Zero-temperature gap `Delta(0)` in Joules, defaulting to weak-coupling `1.7638 k_B T_c`.
    pub fn effective_gap_0(&self) -> f64 {
        self.gap_0.unwrap_or(BCS_GAP_RATIO * K_B * self.tc)
    }

    pub fn read(s: &mut Section) -> Result<Self, MaterialError> {
        let source = read_sources(s, "source")?;
        let tc = s
            .quantity("tc", Dim::Temperature, None)?
            .ok_or_else(|| s.err("tc", "is missing; write tc = \"9.2 K\""))?;
        if !tc.is_finite() || tc <= 0.0 {
            return Err(s.err("tc", format!("critical temperature tc must be positive (found {tc} K)")));
        }
        let raw_gap = s.quantity("gap_0", Dim::Energy, None)?;
        let gap_0 = positive_opt(s, "gap_0", raw_gap, "gap_0")?;
        let raw_ll = s.quantity("lambda_L", Dim::Length, None)?;
        let lambda_l = positive_opt(s, "lambda_L", raw_ll, "lambda_L")?;
        let raw_xi = s.quantity("xi_0", Dim::Length, None)?;
        let xi_0 = positive_opt(s, "xi_0", raw_xi, "xi_0")?;
        let raw_mfp = s.quantity("mean_free_path", Dim::Length, None)?;
        let mean_free_path = positive_opt(s, "mean_free_path", raw_mfp, "mean_free_path")?;
        let raw_sn = s.quantity("sigma_n", Dim::ElectricalConductivity, None)?;
        let sigma_n = positive_opt(s, "sigma_n", raw_sn, "sigma_n")?;
        let raw_dg = s.quantity("dynes_gamma", Dim::Energy, None)?;
        let dynes_gamma = nonneg_opt(s, "dynes_gamma", raw_dg, "dynes_gamma")?;
        let raw_rr = s.quantity("r_res", Dim::SheetResistance, None)?;
        let r_res = nonneg_opt(s, "r_res", raw_rr, "r_res")?;
        let raw_xqp = s.quantity("x_qp_floor", Dim::Dimensionless, None)?;
        let x_qp_floor = nonneg_opt(s, "x_qp_floor", raw_xqp, "x_qp_floor")?;
        let raw_nt = s.quantity("nominal_thickness", Dim::Length, None)?;
        let nominal_thickness = positive_opt(s, "nominal_thickness", raw_nt, "nominal_thickness")?;
        let film_model = match s.string("film_model")? {
            Some(ref text) => Some(FilmModel::parse(text).ok_or_else(|| {
                s.err(
                    "film_model",
                    format!("unknown film_model '{text}'; expected \"kerr\", \"two_sided\" or \"thin_limit\""),
                )
            })?),
            None => None,
        };
        s.finish(SC_KEYS)?;
        Ok(Self {
            source,
            tc,
            gap_0,
            lambda_l,
            xi_0,
            mean_free_path,
            sigma_n,
            dynes_gamma,
            r_res,
            x_qp_floor,
            nominal_thickness,
            film_model,
        })
    }

    pub fn write(&self, tbl: &mut toml::Table) {
        if let Some(v) = write_sources(&self.source) {
            tbl.insert("source".into(), v);
        }
        tbl.insert("tc".into(), format!("{} K", self.tc).into());
        if let Some(g) = self.gap_0 {
            tbl.insert("gap_0".into(), format!("{g} J").into());
        }
        if let Some(ll) = self.lambda_l {
            tbl.insert("lambda_L".into(), format!("{ll} m").into());
        }
        if let Some(xi) = self.xi_0 {
            tbl.insert("xi_0".into(), format!("{xi} m").into());
        }
        if let Some(ell) = self.mean_free_path {
            tbl.insert("mean_free_path".into(), format!("{ell} m").into());
        }
        if let Some(sn) = self.sigma_n {
            tbl.insert("sigma_n".into(), format!("{sn} S/m").into());
        }
        if let Some(dg) = self.dynes_gamma {
            tbl.insert("dynes_gamma".into(), format!("{dg} J").into());
        }
        if let Some(rr) = self.r_res {
            tbl.insert("r_res".into(), format!("{rr} Ohm").into());
        }
        if let Some(xqp) = self.x_qp_floor {
            tbl.insert("x_qp_floor".into(), toml::Value::Float(xqp));
        }
        if let Some(d) = self.nominal_thickness {
            tbl.insert("nominal_thickness".into(), format!("{d} m").into());
        }
        if let Some(fm) = self.film_model {
            tbl.insert("film_model".into(), fm.as_str().into());
        }
    }
}

/// Result of a BCS complex conductivity and thin-film surface impedance evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BcsEval {
    /// Complex conductivity `sigma = sigma_1 - i sigma_2` in S/m (`exp(+i omega t)` convention).
    pub sigma: Complex64,
    /// Normal-state electrical conductivity `sigma_n` in S/m.
    pub sigma_n: f64,
    /// Superconducting energy gap `Delta(T)` at the evaluated temperature in J.
    pub delta_j: f64,
    /// Effective magnetic penetration depth `lambda_eff = (mu_0 omega sigma_2)^(-1/2)` in m.
    pub lambda_eff: f64,
    /// Complex sheet surface impedance `Z_s = R_s + i omega L_k_sq` in Ohm/sq.
    pub z_s: Complex64,
    /// Sheet surface resistance `R_s = Re(Z_s)` in Ohm/sq.
    pub r_s: f64,
    /// Series sheet kinetic inductance `L_k_sq = Im(Z_s) / omega` in H/sq.
    pub l_k_sq: f64,
    /// Equivalent parallel sheet inductance `L_eff = L_k_sq + R_s^2 / (omega^2 L_k_sq)` in H/sq.
    pub l_eff: f64,
    /// Equivalent parallel sheet conductance `G_eff = R_s / (R_s^2 + omega^2 L_k_sq^2)` in S*sq.
    pub g_eff: f64,
}

fn gl32_integrate<F: Fn(f64) -> (f64, f64)>(a: f64, b: f64, f: &F) -> (f64, f64) {
    if b <= a {
        return (0.0, 0.0);
    }
    let mid = 0.5 * (a + b);
    let half = 0.5 * (b - a);
    let (mut s0, mut s1) = (0.0, 0.0);
    for (&x, &w) in GL32_X.iter().zip(GL32_W.iter()) {
        let dx = half * x;
        let (a0, a1) = f(mid - dx);
        let (b0, b1) = f(mid + dx);
        s0 += w * (a0 + b0);
        s1 += w * (a1 + b1);
    }
    (half * s0, half * s1)
}

/// Evaluates `I(alpha) = 2 * int_0^inf du / (exp(alpha * cosh(u)) + 1)` and `dI/dalpha`.
fn bcs_gap_thermal_integral(alpha: f64) -> (f64, f64) {
    let u_max = (1.0 + 42.0 / alpha.max(1e-12)).acosh();
    let panels = [0.0, 0.04 * u_max, 0.12 * u_max, 0.28 * u_max, 0.5 * u_max, 0.75 * u_max, u_max];
    let (mut val, mut dval) = (0.0, 0.0);
    for w in panels.windows(2) {
        let (v, dv) = gl32_integrate(w[0], w[1], &|u| {
            let c = u.cosh();
            let ac = alpha * c;
            if ac > 500.0 {
                return (0.0, 0.0);
            }
            let e = ac.exp();
            let inv = 1.0 / (e + 1.0);
            let f0 = 2.0 * inv;
            let df = -2.0 * c * e * inv * inv;
            (f0, df)
        });
        val += v;
        dval += dv;
    }
    (val, dval)
}

/// Solves the exact self-consistent BCS gap equation for `Delta(T)` in Joules.
///
/// Uses the hyperbolic substitution `E = Delta * cosh(u)` to remove the square-root singularity
/// and solves `ln(Delta(T)/Delta(0)) + 2 int_0^inf du / (exp((pi e^-gamma_E / t) (Delta(T)/Delta(0)) cosh(u)) + 1) = 0`
/// via safeguarded Newton-bisection to `< 1e-13` relative accuracy. Returns `0.0` for `T >= T_c`.
pub fn bcs_gap(tc: f64, gap_0_j: Option<f64>, temperature: f64) -> Result<f64, String> {
    if !tc.is_finite() || tc <= 0.0 {
        return Err(format!("critical temperature tc must be positive and finite (found {tc})"));
    }
    if !temperature.is_finite() || temperature < 0.0 {
        return Err(format!("temperature must be non-negative and finite (found {temperature})"));
    }
    let delta_0 = gap_0_j.unwrap_or(BCS_GAP_RATIO * K_B * tc);
    if !delta_0.is_finite() || delta_0 <= 0.0 {
        return Err(format!("zero-temperature gap gap_0 must be positive and finite (found {delta_0})"));
    }
    if temperature >= tc {
        return Ok(0.0);
    }
    let t = temperature / tc;
    if t <= 0.02 {
        return Ok(delta_0);
    }

    let scale = BCS_EXACT_PI_EXP_GAMMA / t;
    let mut lo = 1e-15_f64;
    let mut hi = 1.0_f64;
    let mut delta = (1.74 * (1.0 / t - 1.0).sqrt()).tanh().clamp(1e-9, 1.0 - 1e-15);

    for _ in 0..40 {
        let alpha = scale * delta;
        let (ith, dith_dalpha) = bcs_gap_thermal_integral(alpha);
        let f = delta.ln() + ith;
        let df = 1.0 / delta + scale * dith_dalpha;
        if f < 0.0 {
            lo = lo.max(delta);
        } else {
            hi = hi.min(delta);
        }
        if f.abs() < 1e-15 || (hi - lo) <= 1e-14 * delta {
            break;
        }
        let next = delta - f / df;
        if df > 0.0 && next > lo && next < hi {
            delta = next;
        } else {
            delta = 0.5 * (lo + hi);
        }
    }
    Ok(delta_0 * delta)
}

#[inline]
fn fermi(e_j: f64, kt_j: f64) -> f64 {
    if kt_j <= 0.0 {
        if e_j > 0.0 {
            0.0
        } else if e_j < 0.0 {
            1.0
        } else {
            0.5
        }
    } else {
        let x = e_j / kt_j;
        if x > 500.0 {
            0.0
        } else if x < -500.0 {
            1.0
        } else {
            1.0 / (x.exp() + 1.0)
        }
    }
}

/// Stable difference `f(e1) - f(e2)` for `e2 = e1 + hw > e1` avoiding cancellation at large `e1 / kT`.
#[inline]
fn fermi_diff(e1: f64, e2: f64, kt_j: f64) -> f64 {
    if kt_j <= 0.0 {
        return if e1 < 0.0 && e2 > 0.0 {
            1.0
        } else if (e1 == 0.0 && e2 > 0.0) || (e1 < 0.0 && e2 == 0.0) {
            0.5
        } else {
            0.0
        };
    }
    let x1 = e1 / kt_j;
    let x2 = e2 / kt_j;
    if x1 > 500.0 || x2 < -500.0 {
        return 0.0;
    }
    let dx = (e2 - e1) / kt_j;
    if x1 > 30.0 {
        let e_x1 = (-x1).exp();
        let e_x2 = (-x2).exp();
        (e_x1 - e_x2) / ((1.0 + e_x1) * (1.0 + e_x2))
    } else if x1 >= 0.0 && dx < 1e-4 {
        let ex1 = x1.exp();
        let edx = dx.exp_m1();
        (ex1 * edx) / ((1.0 + ex1) * (1.0 + ex1 * (1.0 + edx)))
    } else {
        fermi(e1, kt_j) - fermi(e2, kt_j)
    }
}

/// Evaluates `(sigma_1, sigma_2)` in S/m using Zimmermann et al. (1991) when `gamma_s = Some(hbar / tau)`
/// or the exact Mattis-Bardeen dirty limit (`gamma_s = None`).
fn evaluate_bcs_conductivity(
    sigma_n: f64,
    delta_j: f64,
    hw_j: f64,
    kt_j: f64,
    gamma_s: Option<f64>,
) -> (f64, f64) {
    if delta_j <= 0.0 {
        return match gamma_s {
            Some(gs) if gs > 0.0 => {
                let wt = hw_j / gs;
                let denom = 1.0 + wt * wt;
                (sigma_n / denom, sigma_n * wt / denom)
            }
            _ => (sigma_n, 0.0),
        };
    }

    // ----- Interval 1: E_1 in [Delta, inf) via hyperbolic substitution E_1 = Delta * cosh(u) -----
    let (i1_re, i1_im) = {
        let u_max = match gamma_s {
            Some(gs) => (1.0 + 1e9 * (delta_j + gs + hw_j) / delta_j).acosh(),
            None => {
                if kt_j <= 1e-4 * delta_j {
                    0.0
                } else {
                    (1.0 + 45.0 * kt_j / delta_j).acosh()
                }
            }
        };
        if u_max <= 0.0 {
            (0.0, 0.0)
        } else {
            let u_w = (1.0 + hw_j / delta_j).acosh();
            let u_g = gamma_s.map_or(1.0, |gs| (1.0 + gs / delta_j).acosh());
            let u_t = if kt_j > 0.0 { (1.0 + kt_j / delta_j).acosh() } else { 0.1 };
            let mut breaks = vec![
                0.0,
                0.2 * u_w,
                u_w,
                3.0 * u_w,
                0.5 * u_g,
                u_g,
                u_t,
                0.05,
                0.2,
                0.6,
                1.5,
                3.0,
                6.0,
                10.0,
                15.0,
                u_max,
            ];
            breaks.retain(|&u| u >= 0.0 && u <= u_max && u.is_finite());
            breaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
            breaks.dedup_by(|a, b| (*a - *b).abs() < 1e-6);

            let (mut s_re, mut s_im) = (0.0, 0.0);
            for w in breaks.windows(2) {
                let (re, im) = gl32_integrate(w[0], w[1], &|u| {
                    let (sh, ch) = (u.sinh(), u.cosh());
                    let e1 = delta_j * ch;
                    let r1 = delta_j * sh;
                    let e2 = e1 + hw_j;
                    let r2 = ((e2 - delta_j) * (e2 + delta_j)).max(0.0).sqrt();
                    if r2 <= 0.0 {
                        return (0.0, 0.0);
                    }
                    let df = fermi_diff(e1, e2, kt_j);
                    let e1e2_d2 = e1 * e2 + delta_j * delta_j;
                    let r1r2 = r1 * r2;
                    let gp_r1 = (e1e2_d2 + r1r2) / r2;
                    let gm_r1 = (delta_j * delta_j * (e1 + e2) * (e1 + e2)) / (r2 * (e1e2_d2 + r1r2));

                    match gamma_s {
                        None => {
                            let re = 4.0 * df * (e1e2_d2 / r2);
                            (re, 0.0)
                        }
                        Some(gs) => {
                            let gs2 = gs * gs;
                            let r_plus = r1 + r2;
                            let r_minus = hw_j * (e1 + e2) / r_plus;
                            let d_minus = r_minus * r_minus + gs2;
                            let d_plus = r_plus * r_plus + gs2;
                            let re = 2.0 * df * gs2 * (gp_r1 / d_minus + gm_r1 / d_plus);
                            let tanh_sum = 2.0 * (1.0 - fermi(e1, kt_j) - fermi(e2, kt_j));
                            let im = gs * (2.0 * df * gp_r1 * r_minus / d_minus - tanh_sum * gm_r1 * r_plus / d_plus);
                            (re, im)
                        }
                    }
                });
                s_re += re;
                s_im += im;
            }
            (s_re, s_im)
        }
    };

    // ----- Sub-interval 2a: E_1 in [max(-Delta, Delta - hw), Delta] -----
    let e_min_2a = (-delta_j).max(delta_j - hw_j);
    let i2a_im = if delta_j > e_min_2a {
        let mid = 0.5 * (e_min_2a + delta_j);
        let half = 0.5 * (delta_j - e_min_2a);
        let panels = [0.0, 0.25 * PI, 0.5 * PI, 0.75 * PI, PI];
        let mut s_im = 0.0;
        for w in panels.windows(2) {
            let (_, im) = gl32_integrate(w[0], w[1], &|theta| {
                let cos_t = theta.cos();
                let sin_t = theta.sin();
                let e1 = mid + half * cos_t;
                let e2 = e1 + hw_j;
                let d_minus_e1 = (half * (1.0 - cos_t)).max(0.0);
                let e1_minus_emin = (half * (1.0 + cos_t)).max(0.0);
                let s1 = (d_minus_e1 * (delta_j + e1)).max(0.0).sqrt();
                let r2 = ((e2 - delta_j).max(0.0) * (e2 + delta_j)).sqrt();

                let jac_over_s1r2 = if hw_j <= 2.0 * delta_j {
                    1.0 / ((delta_j + e1) * (e2 + delta_j)).sqrt()
                } else {
                    1.0 / r2
                };
                let _ = (sin_t, e1_minus_emin);

                let tanh2 = 1.0 - 2.0 * fermi(e2, kt_j);
                let e1e2_d2 = e1 * e2 + delta_j * delta_j;
                let im = match gamma_s {
                    None => 2.0 * tanh2 * e1e2_d2 * jac_over_s1r2,
                    Some(gs) => {
                        let s1_gs = s1 + gs;
                        let num = s1 * r2 * r2 + e1e2_d2 * s1_gs;
                        let den = r2 * r2 + s1_gs * s1_gs;
                        2.0 * tanh2 * gs * (num / den) * jac_over_s1r2
                    }
                };
                (0.0, im)
            });
            s_im += im;
        }
        s_im
    } else {
        0.0
    };

    // ----- Sub-interval 2b (only when hw > 2 * Delta): E_1 in [Delta - hw, -Delta] -----
    let (i2b_re, i2b_im) = if hw_j > 2.0 * delta_j {
        let e_min = delta_j - hw_j;
        let e_max = -delta_j;
        let mid = 0.5 * (e_min + e_max);
        let half = 0.5 * (e_max - e_min);
        let panels = [0.0, 0.25 * PI, 0.5 * PI, 0.75 * PI, PI];
        let (mut s_re, mut s_im) = (0.0, 0.0);
        for w in panels.windows(2) {
            let (re, im) = gl32_integrate(w[0], w[1], &|theta| {
                let cos_t = theta.cos();
                let e1 = mid + half * cos_t;
                let e2 = e1 + hw_j;
                let emax_minus_e1 = (half * (1.0 - cos_t)).max(0.0);
                let e1_minus_emin = (half * (1.0 + cos_t)).max(0.0);
                let r1 = (emax_minus_e1 * (-e1 + delta_j)).sqrt();
                let r2 = (e1_minus_emin * (e2 + delta_j)).sqrt();
                let jac_over_r1r2 = 1.0 / ((-e1 + delta_j) * (e2 + delta_j)).sqrt();
                let tanh2 = 1.0 - 2.0 * fermi(e2, kt_j);
                let e1e2_d2 = e1 * e2 + delta_j * delta_j;
                let r1r2 = r1 * r2;
                match gamma_s {
                    None => {
                        let re = 2.0 * tanh2 * (-e1e2_d2) * jac_over_r1r2;
                        (re, 0.0)
                    }
                    Some(gs) => {
                        let gs2 = gs * gs;
                        let r_plus = r1 + r2;
                        let r_minus = r2 - r1;
                        let d_plus = r_plus * r_plus + gs2;
                        let d_minus = r_minus * r_minus + gs2;
                        let one_minus_g = r1r2 - e1e2_d2;
                        let minus_one_minus_g = -(r1r2 + e1e2_d2);
                        let re = tanh2 * gs2 * (one_minus_g / d_plus + minus_one_minus_g / d_minus) * jac_over_r1r2;
                        let im = tanh2 * gs * (one_minus_g * r_plus / d_plus + minus_one_minus_g * r_minus / d_minus) * jac_over_r1r2;
                        (re, im)
                    }
                }
            });
            s_re += re;
            s_im += im;
        }
        (s_re, s_im)
    } else {
        (0.0, 0.0)
    };

    let scale = sigma_n / (2.0 * hw_j);
    let sigma_1 = (scale * (i1_re + i2b_re)).max(0.0);
    let sigma_2 = (scale * (i1_im + i2a_im + i2b_im)).max(0.0);
    (sigma_1, sigma_2)
}

/// Dynes sub-gap pair-breaking real conductivity `sigma_1,Dynes >= 0` integrated strictly over the
/// sub-gap channel (`|E_1| < Delta` or `|E_2| < Delta`) to avoid any double counting with the
/// Zimmermann/Mattis-Bardeen coherence integral.
fn dynes_subgap_sigma_1(sigma_n: f64, delta_j: f64, hw_j: f64, kt_j: f64, gamma_d: f64) -> f64 {
    if delta_j <= 0.0 || gamma_d <= 0.0 {
        return 0.0;
    }
    let coeffs = |e: f64| -> (f64, f64) {
        let z = Complex64::new(e, -gamma_d);
        let mut s = (z * z - delta_j * delta_j).sqrt();
        if (z / s).re < 0.0 {
            s = -s;
        }
        let gr = (z / s).re.max(0.0);
        let ga = (Complex64::new(delta_j, 0.0) / s).re;
        (gr, ga)
    };
    let e_lo = -delta_j - hw_j;
    let e_hi = delta_j;
    let mut pts = vec![e_lo, -delta_j, (-hw_j).clamp(e_lo, e_hi), 0.0, (delta_j - hw_j).clamp(e_lo, e_hi), e_hi];
    pts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    pts.dedup_by(|a, b| (*a - *b).abs() < 1e-35);
    let mut integral = 0.0;
    for w in pts.windows(2) {
        let (re, _) = gl32_integrate(w[0], w[1], &|e1| {
            let e2 = e1 + hw_j;
            if e1.abs() >= delta_j && e2.abs() >= delta_j {
                return (0.0, 0.0);
            }
            let df = fermi_diff(e1, e2, kt_j);
            let (gr1, ga1) = coeffs(e1);
            let (gr2, ga2) = coeffs(e2);
            ((df * (gr1 * gr2 + ga1 * ga2)).max(0.0), 0.0)
        });
        integral += re;
    }
    (sigma_n / hw_j) * integral.max(0.0)
}

/// Effective quasiparticle temperature `T_eff` corresponding to thermal normalized quasiparticle
/// density `sqrt(2 pi k_B T / Delta) * exp(-Delta / (k_B T)) = x_qp`.
fn effective_temp_for_xqp(delta_j: f64, x_qp: f64) -> f64 {
    if delta_j <= 0.0 || x_qp <= 0.0 {
        return 0.0;
    }
    let target = x_qp.min(0.4);
    let u = (target / (2.0 * PI).sqrt()).ln();
    // Solve ln(sqrt(t)) - 1/t = ln(target / sqrt(2 pi)), with t = k_B T / Delta in (0, 0.5]
    let mut lo = 1e-4_f64;
    let mut hi = 0.6_f64;
    let mut t = (1.0 / (-u).max(2.0)).clamp(lo, hi);
    for _ in 0..40 {
        let g = 0.5 * t.ln() - 1.0 / t - u;
        if g < 0.0 {
            lo = lo.max(t);
        } else {
            hi = hi.min(t);
        }
        let dg = 0.5 / t + 1.0 / (t * t);
        let next = t - g / dg;
        t = if next > lo && next < hi { next } else { 0.5 * (lo + hi) };
        if (hi - lo) < 1e-12 * t {
            break;
        }
    }
    t * delta_j / K_B
}

/// Evaluates arbitrary-purity BCS complex conductivity and thin-film surface impedance from explicit parameters.
#[allow(clippy::too_many_arguments)]
pub fn surface_impedance_from_params(
    tc: f64,
    gap_0_j: Option<f64>,
    lambda_l: Option<f64>,
    xi_0: Option<f64>,
    mean_free_path: Option<f64>,
    sigma_n: Option<f64>,
    freq_hz: f64,
    temperature: f64,
    thickness: f64,
    film_model: FilmModel,
) -> Result<BcsEval, String> {
    surface_impedance_full(
        tc,
        gap_0_j,
        lambda_l,
        xi_0,
        mean_free_path,
        sigma_n,
        None,
        None,
        None,
        freq_hz,
        temperature,
        thickness,
        film_model,
    )
}

#[allow(clippy::too_many_arguments)]
fn surface_impedance_full(
    tc: f64,
    gap_0_j: Option<f64>,
    lambda_l: Option<f64>,
    xi_0: Option<f64>,
    mean_free_path: Option<f64>,
    sigma_n: Option<f64>,
    dynes_gamma: Option<f64>,
    r_res: Option<f64>,
    x_qp_floor: Option<f64>,
    freq_hz: f64,
    temperature: f64,
    thickness: f64,
    film_model: FilmModel,
) -> Result<BcsEval, String> {
    if !freq_hz.is_finite() || freq_hz <= 0.0 {
        return Err(format!("frequency freq_hz must be positive and finite (found {freq_hz})"));
    }
    if !thickness.is_finite() || thickness <= 0.0 {
        return Err(format!("film thickness must be positive and finite (found {thickness})"));
    }
    let delta_0 = gap_0_j.unwrap_or(BCS_GAP_RATIO * K_B * tc);
    let delta_j = bcs_gap(tc, Some(delta_0), temperature)?;

    for (name, val) in [("lambda_L", lambda_l), ("xi_0", xi_0), ("mean_free_path", mean_free_path), ("sigma_n", sigma_n)] {
        if let Some(v) = val {
            if !v.is_finite() || v <= 0.0 {
                return Err(format!("{name} must be positive and finite (found {v})"));
            }
        }
    }

    let (sigma_n_val, gamma_s) = match (sigma_n, lambda_l, xi_0, mean_free_path) {
        (Some(sn), Some(ll), Some(xi), mfp) => {
            let v_f = PI * delta_0 * xi / HBAR;
            let ell = mfp.unwrap_or(sn * MU0 * ll * ll * v_f);
            (sn, Some(PI * delta_0 * xi / ell))
        }
        (None, Some(ll), Some(xi), Some(ell)) => {
            let v_f = PI * delta_0 * xi / HBAR;
            let sn = ell / (MU0 * ll * ll * v_f);
            (sn, Some(PI * delta_0 * xi / ell))
        }
        (None, Some(ll), xi_opt, None) => {
            let xi = xi_opt.unwrap_or(100e-9);
            let v_f = PI * delta_0 * xi / HBAR;
            let ell = 1e6 * xi;
            let sn = ell / (MU0 * ll * ll * v_f);
            (sn, Some(PI * delta_0 * xi / ell))
        }
        (Some(sn), None, Some(xi), Some(ell)) => (sn, Some(PI * delta_0 * xi / ell)),
        (Some(sn), _, _, _) => (sn, None),
        _ => {
            return Err(
                "superconducting surface impedance requires sigma_n or (lambda_L and xi_0 with sigma_n or mean_free_path)".into(),
            );
        }
    };

    let omega = 2.0 * PI * freq_hz;
    let hw_j = HBAR * omega;
    let kt_j = K_B * temperature;

    let (mut sigma_1, sigma_2) = evaluate_bcs_conductivity(sigma_n_val, delta_j, hw_j, kt_j, gamma_s);

    // Non-equilibrium quasiparticle density floor adds dissipative sigma_1 evaluated at T_qp(x_qp_floor).
    if let Some(xqp) = x_qp_floor.filter(|&x| x > 0.0 && delta_j > 0.0) {
        let t_eff = effective_temp_for_xqp(delta_j, xqp);
        if t_eff > temperature {
            let (s1_floor, _) = evaluate_bcs_conductivity(sigma_n_val, delta_j, hw_j, K_B * t_eff, gamma_s);
            sigma_1 = sigma_1.max(s1_floor);
        }
    }

    // Dynes sub-gap pair-breaking density of states adds finite sub-gap absorption.
    if let Some(gd) = dynes_gamma.filter(|&g| g > 0.0 && delta_j > 0.0) {
        sigma_1 += dynes_subgap_sigma_1(sigma_n_val, delta_j, hw_j, kt_j, gd);
    }

    let sigma = Complex64::new(sigma_1, -sigma_2);

    let lambda_eff = if sigma_2 > 0.0 {
        1.0 / (MU0 * omega * sigma_2).sqrt()
    } else {
        f64::INFINITY
    };

    let i_omega_mu0 = Complex64::new(0.0, omega * MU0);
    let gamma = (i_omega_mu0 * sigma).sqrt();

    let coth_scaled = |z: Complex64| -> Complex64 {
        let z = if z.re < 0.0 { -z } else { z };
        if z.norm() < 1e-4 {
            let z2 = z * z;
            Complex64::new(1.0, 0.0) + z2 / 3.0 - (z2 * z2) / 45.0
        } else if z.re > 300.0 {
            z
        } else {
            let e2 = (-2.0 * z).exp();
            z * (Complex64::new(1.0, 0.0) + e2) / (Complex64::new(1.0, 0.0) - e2)
        }
    };

    let inv_sigma_d = Complex64::new(1.0, 0.0) / (sigma * thickness);
    let mut z_s = match film_model {
        FilmModel::ThinLimit => inv_sigma_d,
        FilmModel::Kerr => inv_sigma_d * coth_scaled(gamma * thickness),
        FilmModel::TwoSided => inv_sigma_d * coth_scaled(0.5 * gamma * thickness),
    };

    let r_s = z_s.re.max(0.0) + r_res.unwrap_or(0.0).max(0.0);
    z_s = Complex64::new(r_s, z_s.im);
    let l_k_sq = z_s.im / omega;
    let l_eff = if l_k_sq.abs() > 0.0 {
        l_k_sq + (r_s * r_s) / (omega * omega * l_k_sq)
    } else {
        0.0
    };
    let denom = r_s * r_s + (omega * l_k_sq) * (omega * l_k_sq);
    let g_eff = if denom > 0.0 { r_s / denom } else { 0.0 };

    Ok(BcsEval {
        sigma,
        sigma_n: sigma_n_val,
        delta_j,
        lambda_eff,
        z_s,
        r_s,
        l_k_sq,
        l_eff,
        g_eff,
    })
}

/// Evaluates arbitrary-purity BCS complex conductivity and thin-film surface impedance for a `[superconducting]` block,
/// including sub-gap `dynes_gamma`, `r_res`, and `x_qp_floor` loss channels when present.
pub fn bcs_surface_impedance(
    sc: &Superconducting,
    freq_hz: f64,
    temperature: f64,
    thickness: f64,
    film_model: FilmModel,
) -> Result<BcsEval, String> {
    surface_impedance_full(
        sc.tc,
        sc.gap_0,
        sc.lambda_l,
        sc.xi_0,
        sc.mean_free_path,
        sc.sigma_n,
        sc.dynes_gamma,
        sc.r_res,
        sc.x_qp_floor,
        freq_hz,
        temperature,
        thickness,
        film_model,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bcs_gap_self_consistency_and_bounds() {
        let tc = 9.2;
        let d0 = BCS_GAP_RATIO * K_B * tc;
        assert_eq!(bcs_gap(tc, None, 0.0).unwrap(), d0);
        assert_eq!(bcs_gap(tc, None, tc).unwrap(), 0.0);
        assert_eq!(bcs_gap(tc, None, 12.0).unwrap(), 0.0);
        for k in 1..20 {
            let t = (k as f64 * 0.05) * tc;
            let d = bcs_gap(tc, None, t).unwrap();
            assert!(d > 0.0 && d <= d0);
            let delta_norm = d / d0;
            let alpha = BCS_EXACT_PI_EXP_GAMMA * delta_norm / (t / tc);
            let (ith, _) = bcs_gap_thermal_integral(alpha);
            assert!(
                (delta_norm.ln() + ith).abs() < 1e-13,
                "self-consistency residual at T={t}: {}",
                delta_norm.ln() + ith
            );
        }
    }

    #[test]
    fn dirty_and_clean_limits_and_park_nb() {
        let tc = 9.2;
        let gap_0 = 1.395e-3 * EV_TO_J;
        let freq = 7.0e9;
        let omega = 2.0 * PI * freq;

        // Dirty limit: sigma_2(0) / sigma_n -> pi * Delta_0 / (hbar * omega)
        let dirty = surface_impedance_from_params(
            tc,
            Some(gap_0),
            None,
            None,
            None,
            Some(5.0e5),
            1.0e9,
            0.01,
            50e-9,
            FilmModel::ThinLimit,
        )
        .unwrap();
        let want_dirty = 5.0e5 * PI * gap_0 / (HBAR * 2.0 * PI * 1.0e9);
        assert!((-dirty.sigma.im - want_dirty).abs() < 2e-4 * want_dirty);

        // Clean limit: lambda_eff(0) -> lambda_L
        let lambda_l = 33.3e-9;
        let clean = surface_impedance_from_params(
            tc,
            Some(gap_0),
            Some(lambda_l),
            Some(39.0e-9),
            Some(1e5 * 39.0e-9),
            None,
            freq,
            0.01,
            100e-9,
            FilmModel::Kerr,
        )
        .unwrap();
        assert!((clean.lambda_eff - lambda_l).abs() < 1e-4 * lambda_l);

        // Park et al. 2026 Nb (d = 100 nm, f = 7 GHz, T = 10 mK)
        let nb = Superconducting {
            tc,
            gap_0: Some(gap_0),
            lambda_l: Some(lambda_l),
            xi_0: Some(39.0e-9),
            sigma_n: Some(5.5e7),
            ..Default::default()
        };
        let kerr = bcs_surface_impedance(&nb, freq, 0.01, 100e-9, FilmModel::Kerr).unwrap();
        let two = bcs_surface_impedance(&nb, freq, 0.01, 100e-9, FilmModel::TwoSided).unwrap();
        let thin = bcs_surface_impedance(&nb, freq, 0.01, 100e-9, FilmModel::ThinLimit).unwrap();
        assert!(kerr.l_k_sq > two.l_k_sq && two.l_k_sq > thin.l_k_sq && thin.l_k_sq > 0.0);
        assert!((kerr.z_s.im / omega - kerr.l_k_sq).abs() < 1e-25);
    }

    #[test]
    fn subgap_mk_loss_channels_give_finite_positive_resistance() {
        let base = Superconducting {
            tc: 9.2,
            gap_0: Some(1.395e-3 * EV_TO_J),
            lambda_l: Some(33.3e-9),
            xi_0: Some(39.0e-9),
            sigma_n: Some(5.5e7),
            ..Default::default()
        };
        let ev0 = bcs_surface_impedance(&base, 5.0e9, 0.01, 100e-9, FilmModel::Kerr).unwrap();

        let with_rres = Superconducting {
            r_res: Some(12e-9),
            ..base.clone()
        };
        let ev_r = bcs_surface_impedance(&with_rres, 5.0e9, 0.01, 100e-9, FilmModel::Kerr).unwrap();
        assert!(ev_r.r_s >= ev0.r_s + 11.9e-9 && ev_r.r_s.is_finite());

        let with_dynes = Superconducting {
            dynes_gamma: Some(5e-6 * EV_TO_J),
            ..base.clone()
        };
        let ev_d = bcs_surface_impedance(&with_dynes, 5.0e9, 0.01, 100e-9, FilmModel::Kerr).unwrap();
        assert!(ev_d.sigma.re > ev0.sigma.re && ev_d.r_s > ev0.r_s);

        let with_xqp = Superconducting {
            x_qp_floor: Some(1e-5),
            ..base
        };
        let ev_x = bcs_surface_impedance(&with_xqp, 5.0e9, 0.01, 100e-9, FilmModel::Kerr).unwrap();
        assert!(ev_x.sigma.re > ev0.sigma.re && ev_x.r_s > ev0.r_s);
    }
}
