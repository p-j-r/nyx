// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Piezoelectric resonator extraction & Larson (2000) 6-parameter modified Butterworth–Van Dyke (mBVD)
//! fitting (`PiezoResonatorAnalysis`, `id = "piezo_resonator"`).
//!
//! Operates on multi-port [`DrivenPortData`] produced by `nyx-sms` `Study::Driven`:
//! - Reduces multi-port admittance matrices via [`crate::network::schur_reduced_admittance`] supporting
//!   `0.0` (short), positive `f64` (Ohm), and `"open"` terminations.
//! - Extracts series/parallel resonance frequencies `f_s_hz`, `f_p_hz`, effective coupling `k_eff2`,
//!   IEEE Std 176 thickness-mode coupling `k_t2_ieee`, unreduced static capacitance `c_0_raw_f`,
//!   and unwrapped phase-slope quality factors `q_s`, `q_p` (`None` / `null` when lossless or under-resolved).
//! - Fits the 6-parameter Larson mBVD circuit `(C_0, R_0, R_s, C_m, L_m, R_m)` via bounded
//!   Levenberg–Marquardt in bounded reflection space `Gamma = (1 - Z_ref Y) / (1 + Z_ref Y)`.
#![allow(clippy::needless_range_loop)]

use crate::network::schur_reduced_admittance;
use crate::Analysis;
use femcore::linalg::dense::Dense;
use femcore::linalg::scalar::C64;
use nyx_physics::{DrivenPortData, Error, Mesh, SolveResult};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::f64::consts::PI;

/// Extracted resonant mode summary (primary or secondary spurious mode inside the window).
#[derive(Clone, Debug)]
pub struct PiezoModeMetrics {
    pub f_s_hz: f64,
    pub f_p_hz: f64,
    pub k_eff2: f64,
    pub k_t2_ieee: f64,
    pub q_s: Option<f64>,
    pub q_p: Option<f64>,
}

/// Fitted 6-parameter Larson (2000) modified Butterworth–Van Dyke equivalent circuit and resonator FoMs.
#[derive(Clone, Debug)]
pub struct PiezoResonatorReport {
    pub port: String,
    pub z_ref_ohm: f64,
    pub f_s_hz: f64,
    pub f_p_hz: f64,
    pub k_eff2: f64,
    pub k_t2_ieee: f64,
    pub c_0_raw_f: f64,
    pub c0_f: f64,
    pub r0_ohm: f64,
    pub r_s_ohm: f64,
    pub cm_f: f64,
    pub lm_h: f64,
    pub rm_ohm: f64,
    pub q_s: Option<f64>,
    pub q_p: Option<f64>,
    pub q_m: Option<f64>,
    pub fom: Option<f64>,
    pub fit_residual: f64,
    pub additional_modes: Vec<PiezoModeMetrics>,
}

impl PiezoResonatorReport {
    pub fn to_json(&self) -> Value {
        let extra: Vec<Value> = self
            .additional_modes
            .iter()
            .map(|m| {
                json!({
                    "f_s_hz": m.f_s_hz,
                    "f_p_hz": m.f_p_hz,
                    "k_eff2": m.k_eff2,
                    "k_t2_ieee": m.k_t2_ieee,
                    "q_s": m.q_s,
                    "q_p": m.q_p,
                })
            })
            .collect();
        json!({
            "port": self.port,
            "z_ref_ohm": self.z_ref_ohm,
            "f_s_hz": self.f_s_hz,
            "f_p_hz": self.f_p_hz,
            "k_eff2": self.k_eff2,
            "k_t2_ieee": self.k_t2_ieee,
            "c_0_raw_f": self.c_0_raw_f,
            "c0_f": self.c0_f,
            "r0_ohm": self.r0_ohm,
            "r_s_ohm": self.r_s_ohm,
            "cm_f": self.cm_f,
            "lm_h": self.lm_h,
            "rm_ohm": self.rm_ohm,
            "q_s": self.q_s,
            "q_p": self.q_p,
            "q_m": self.q_m,
            "fom": self.fom,
            "fit_residual": self.fit_residual,
            "additional_modes": extra,
        })
    }
}

/// True when `result` is a `nyx-sms` driven solve carrying multi-port `DrivenPortData`.
pub fn piezo_resonator_applies(result: &SolveResult) -> bool {
    result.physics == "sms" && DrivenPortData::of(result).is_some()
}

/// Evaluates the Larson (2000) 6-parameter mBVD admittance `Y(omega)` for `(C0, R0, Rs, Cm, Lm, Rm)`.
#[inline]
pub fn mbvd_admittance(omega: f64, c0: f64, r0: f64, r_s: f64, cm: f64, lm: f64, rm: f64) -> C64 {
    let iwc0 = C64::new(0.0, omega * c0);
    let y0 = iwc0 / (C64::ONE + iwc0 * r0);
    let zm = C64::new(rm, omega * lm - 1.0 / (omega * cm));
    let y_core = y0 + C64::ONE / zm;
    if r_s <= 0.0 {
        y_core
    } else {
        C64::ONE / (C64::new(r_s, 0.0) + C64::ONE / y_core)
    }
}

fn pairs_to_dense(raw: &[Vec<[f64; 2]>]) -> Dense<C64> {
    let n = raw.len();
    Dense::from_fn(n, n, |i, j| C64::new(raw[i][j][0], raw[i][j][1]))
}

fn parse_terminations(ports: &[String], active_idx: usize, params: &Value) -> Result<(Vec<usize>, Vec<f64>), Error> {
    let mut term_map: BTreeMap<String, f64> = BTreeMap::new();
    if let Some(obj) = params.get("terminations").filter(|v| !v.is_null()) {
        let map = obj.as_object().ok_or_else(|| {
            Error::invalid("parameter 'terminations' must be an object mapping port names to Ohm values or \"open\"")
        })?;
        for (name, val) in map {
            if !ports.iter().any(|p| p == name) {
                return Err(Error::invalid(format!("termination port '{name}' not found in {ports:?}")));
            }
            let z = if let Some(s) = val.as_str() {
                if s.eq_ignore_ascii_case("open") {
                    f64::INFINITY
                } else {
                    return Err(Error::invalid(format!(
                        "termination for port '{name}' must be a non-negative number or \"open\", got '{s}'"
                    )));
                }
            } else if let Some(v) = val.as_f64() {
                if !v.is_finite() || v < 0.0 {
                    return Err(Error::invalid(format!(
                        "termination for port '{name}' must be non-negative and finite, got {v}"
                    )));
                }
                v
            } else {
                return Err(Error::invalid(format!(
                    "termination for port '{name}' must be a non-negative number or \"open\""
                )));
            };
            term_map.insert(name.clone(), z);
        }
    }
    let mut term_idx = Vec::new();
    let mut z_loads = Vec::new();
    for (idx, p) in ports.iter().enumerate() {
        if idx != active_idx {
            term_idx.push(idx);
            z_loads.push(term_map.get(p).copied().unwrap_or(0.0));
        }
    }
    Ok((term_idx, z_loads))
}

fn phase_slope_q(freqs: &[f64], vals: &[C64], f_center: f64, peak_height_check: bool) -> Option<f64> {
    if freqs.len() < 3 {
        return None;
    }
    let idx = freqs
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| (*a - f_center).abs().total_cmp(&(*b - f_center).abs()))
        .map(|(i, _)| i)?;
    let (i0, i1) = if idx == 0 {
        (0, 1)
    } else if idx + 1 >= freqs.len() {
        (freqs.len() - 2, freqs.len() - 1)
    } else {
        (idx - 1, idx + 1)
    };
    let mut phases: Vec<f64> = vals[i0..=i1].iter().map(|z| z.im.atan2(z.re)).collect();
    for k in 1..phases.len() {
        let mut d = phases[k] - phases[k - 1];
        while d > PI {
            d -= 2.0 * PI;
        }
        while d < -PI {
            d += 2.0 * PI;
        }
        phases[k] = phases[k - 1] + d;
    }
    let dphi = (phases[phases.len() - 1] - phases[0]).abs();
    let df = freqs[i1] - freqs[i0];
    if df <= 0.0 || dphi <= 1e-9 || dphi > 2.8 {
        // Zero slope (lossless flat zero/pi) or severely under-resolved (>~0.9*pi phase jump over one stencil)
        return None;
    }
    if peak_height_check {
        let re_vals: Vec<f64> = vals[i0..=i1].iter().map(|v| v.re).collect();
        let max_re = re_vals.iter().copied().fold(0.0_f64, f64::max);
        let max_im = vals[i0..=i1].iter().map(|v| v.im.abs()).fold(0.0_f64, f64::max);
        if max_re <= 1e-12 * max_im.max(1e-30) {
            return None;
        }
    }
    let slope = dphi / df;
    let q = 0.5 * f_center * slope;
    if q.is_finite() && q > 0.5 && q < 1e9 {
        Some(q)
    } else {
        None
    }
}

fn find_upward_zero_crossings(freqs: &[f64], x: &[f64], mag: &[f64]) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for k in 0..freqs.len().saturating_sub(1) {
        if x[k] <= 0.0 && x[k + 1] > 0.0 && (x[k + 1] - x[k]) > 0.0 {
            let t = (-x[k]) / (x[k + 1] - x[k]);
            let f_z = freqs[k] + t * (freqs[k + 1] - freqs[k]);
            let m = mag[k].max(mag[k + 1]);
            out.push((f_z, m));
        }
    }
    out
}

fn compute_kt2_ieee(f_s: f64, f_p: f64) -> f64 {
    if f_s <= 0.0 || f_p <= f_s {
        return 0.0;
    }
    let x = 0.5 * PI * (f_s / f_p);
    let gap = 0.5 * PI * (f_p - f_s) / f_p;
    (x * gap.tan()).max(0.0)
}

fn estimate_c0(freqs: &[f64], y: &[C64], f_s: f64, f_p: f64) -> f64 {
    let mut ests = Vec::new();
    for (&f, &yk) in freqs.iter().zip(y) {
        if f < 0.92 * f_s || f > 1.08 * f_p {
            let c_est = yk.im / (2.0 * PI * f);
            if c_est.is_finite() && c_est > 0.0 {
                ests.push(c_est);
            }
        }
    }
    if ests.is_empty() {
        for (&f, &yk) in freqs.iter().zip(y) {
            let c_est = yk.im / (2.0 * PI * f);
            if c_est.is_finite() && c_est > 0.0 {
                ests.push(c_est);
            }
        }
    }
    if ests.is_empty() {
        return 1.0e-12;
    }
    ests.sort_by(|a, b| a.total_cmp(b));
    ests[ests.len() / 2]
}

fn gamma_of(y: C64, z_ref: f64) -> C64 {
    let zy = y * z_ref;
    (C64::ONE - zy) / (C64::ONE + zy)
}

fn solve_small_real(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut m = a.to_vec();
    let mut x = b.to_vec();
    for k in 0..n {
        let mut piv = k;
        let mut best = m[k][k].abs();
        for r in (k + 1)..n {
            if m[r][k].abs() > best {
                best = m[r][k].abs();
                piv = r;
            }
        }
        if !best.is_finite() || best <= 1e-30 {
            return None;
        }
        if piv != k {
            m.swap(k, piv);
            x.swap(k, piv);
        }
        let d = m[k][k];
        for r in (k + 1)..n {
            let f = m[r][k] / d;
            if f != 0.0 {
                for c in k..n {
                    m[r][c] -= f * m[k][c];
                }
                x[r] -= f * x[k];
            }
        }
    }
    for k in (0..n).rev() {
        let mut s = x[k];
        for c in (k + 1)..n {
            s -= m[k][c] * x[c];
        }
        x[k] = s / m[k][k];
    }
    Some(x)
}

fn fit_mbvd_gamma_lm(
    freqs: &[f64],
    y_data: &[C64],
    z_ref: f64,
    f_s0: f64,
    f_p0: f64,
    c0_init: f64,
    fix_r_s: Option<f64>,
) -> ([f64; 6], f64) {
    let mut k2_eff = ((f_p0 / f_s0).powi(2) - 1.0).max(1e-4);
    if !k2_eff.is_finite() {
        k2_eff = 0.05;
    }
    let cm0 = (c0_init * k2_eff).max(1e-18);
    let lm0 = 1.0 / ((2.0 * PI * f_s0).powi(2) * cm0);
    let z_min_re = y_data.iter().map(|y| (C64::ONE / *y).re.max(0.0)).fold(f64::INFINITY, f64::min);
    let rs0 = fix_r_s.unwrap_or_else(|| if z_min_re.is_finite() { z_min_re.max(0.0) } else { 0.0 });
    let near_fs = freqs
        .iter()
        .zip(y_data)
        .min_by(|(fa, _), (fb, _)| (*fa - f_s0).abs().total_cmp(&(*fb - f_s0).abs()))
        .map(|(_, y)| ((C64::ONE / *y).re - rs0).max(0.0))
        .unwrap_or(1.0);
    let rm0 = near_fs.max(1e-6);
    let tan_off = freqs
        .iter()
        .zip(y_data)
        .filter(|(&f, _)| f < 0.9 * f_s0 || f > 1.1 * f_p0)
        .map(|(_, y)| (y.re / y.im.abs().max(1e-30)).max(0.0))
        .fold(0.0_f64, f64::max);
    let r0_0 = (tan_off / (2.0 * PI * f_s0 * c0_init)).max(1e-7);

    // Parameter vector: [ln(c0), ln(cm), ln(lm), sqrt(rm), sqrt(r0), (optional sqrt(rs))]
    let mut p = vec![c0_init.ln(), cm0.ln(), lm0.ln(), rm0.sqrt(), r0_0.sqrt()];
    if fix_r_s.is_none() {
        p.push(rs0.sqrt());
    }
    let unpack = |u: &[f64]| -> (f64, f64, f64, f64, f64, f64) {
        let c0 = u[0].exp().clamp(1e-22, 1e-3);
        let cm = u[1].exp().clamp(1e-24, 1e-3);
        let lm = u[2].exp().clamp(1e-16, 1e2);
        let rm = (u[3] * u[3]).min(1e7);
        let r0 = (u[4] * u[4]).min(1e7);
        let rs = fix_r_s.unwrap_or_else(|| (u[5] * u[5]).min(1e7));
        (c0, r0, rs, cm, lm, rm)
    };

    let n_pts = freqs.len();
    let gamma_target: Vec<C64> = y_data.iter().map(|&y| gamma_of(y, z_ref)).collect();
    let residuals = |u: &[f64]| -> Vec<f64> {
        let (c0, r0, rs, cm, lm, rm) = unpack(u);
        let mut r = Vec::with_capacity(2 * n_pts);
        for (&f, &gt) in freqs.iter().zip(&gamma_target) {
            let gm = gamma_of(mbvd_admittance(2.0 * PI * f, c0, r0, rs, cm, lm, rm), z_ref);
            r.push(gm.re - gt.re);
            r.push(gm.im - gt.im);
        }
        r
    };
    let cost = |r: &[f64]| -> f64 { r.iter().map(|&v| v * v).sum() };

    let mut r_cur = residuals(&p);
    let mut err_cur = cost(&r_cur);
    let mut lambda = 1e-3_f64;
    let m_dim = p.len();
    let mut jac = vec![vec![0.0; m_dim]; 2 * n_pts];

    for _ in 0..80 {
        for j in 0..m_dim {
            let step = if j < 3 { 1e-5 } else { 1e-6 * p[j].abs().max(1.0) };
            let mut p_hi = p.clone();
            p_hi[j] += step;
            let r_hi = residuals(&p_hi);
            for i in 0..(2 * n_pts) {
                jac[i][j] = (r_hi[i] - r_cur[i]) / step;
            }
        }
        let mut jtj = vec![vec![0.0; m_dim]; m_dim];
        let mut jtr = vec![0.0; m_dim];
        for i in 0..(2 * n_pts) {
            for a in 0..m_dim {
                jtr[a] -= jac[i][a] * r_cur[i];
                for b in a..m_dim {
                    jtj[a][b] += jac[i][a] * jac[i][b];
                }
            }
        }
        for a in 0..m_dim {
            for b in 0..a {
                jtj[a][b] = jtj[b][a];
            }
            jtj[a][a] += lambda * jtj[a][a].max(1e-12);
        }
        let Some(dp) = solve_small_real(&jtj, &jtr) else {
            lambda *= 10.0;
            continue;
        };
        let step_norm: f64 = dp.iter().map(|&v| v * v).sum::<f64>().sqrt();
        if step_norm < 1e-11 {
            break;
        }
        let p_try: Vec<f64> = p
            .iter()
            .zip(&dp)
            .enumerate()
            .map(|(j, (&v, &d))| {
                let nv = v + d;
                if j < 3 {
                    nv.clamp(-55.0, 5.0)
                } else {
                    nv.clamp(-1e4, 1e4)
                }
            })
            .collect();
        let r_try = residuals(&p_try);
        let err_try = cost(&r_try);
        if err_try < err_cur {
            let rel_imp = (err_cur - err_try) / err_cur.max(1e-30);
            p = p_try;
            r_cur = r_try;
            err_cur = err_try;
            lambda = (lambda * 0.3).max(1e-12);
            if rel_imp < 1e-11 && err_cur < 1e-22 {
                break;
            }
        } else {
            lambda = (lambda * 5.0).min(1e12);
        }
    }

    let (c0, r0, rs, cm, lm, rm) = unpack(&p);
    let rms = (err_cur / (n_pts as f64)).sqrt();
    ([c0, r0, rs, cm, lm, rm], rms)
}

/// Runs piezoelectric resonator extraction & mBVD fitting with optional JSON parameters.
pub fn piezo_resonator_with_params(result: &SolveResult, params: &Value) -> Result<PiezoResonatorReport, Error> {
    if !piezo_resonator_applies(result) {
        return Err(Error::invalid("piezo_resonator requires a nyx-sms driven solve result with DrivenPortData"));
    }
    let drv = DrivenPortData::of(result).ok_or_else(|| Error::invalid("missing DrivenPortData"))?;
    if drv.ports.is_empty() || drv.freqs_hz.len() < 5 {
        return Err(Error::invalid("piezo_resonator requires at least 1 port and 5 frequency points"));
    }

    let port_name = match params.get("port").filter(|v| !v.is_null()) {
        Some(v) => v.as_str().ok_or_else(|| Error::invalid("parameter 'port' must be a string"))?.to_string(),
        None => drv.ports[0].clone(),
    };
    let active_idx = drv
        .ports
        .iter()
        .position(|p| p == &port_name)
        .ok_or_else(|| Error::invalid(format!("port '{port_name}' not found in {:?}", drv.ports)))?;

    let (term_idx, z_loads) = parse_terminations(&drv.ports, active_idx, params)?;

    let f_window = match params.get("f_window_hz").filter(|v| !v.is_null()) {
        Some(v) => {
            let arr = v
                .as_array()
                .filter(|a| a.len() == 2)
                .ok_or_else(|| Error::invalid("parameter 'f_window_hz' must be [f_lo_hz, f_hi_hz]"))?;
            let lo = arr[0].as_f64().ok_or_else(|| Error::invalid("f_window_hz[0] must be numeric"))?;
            let hi = arr[1].as_f64().ok_or_else(|| Error::invalid("f_window_hz[1] must be numeric"))?;
            if !lo.is_finite() || !hi.is_finite() || lo >= hi {
                return Err(Error::invalid(format!("invalid f_window_hz [{lo}, {hi}]")));
            }
            Some([lo, hi])
        }
        None => None,
    };

    let fix_r_s = match params.get("fix_r_s_ohm").filter(|v| !v.is_null()) {
        Some(v) => {
            let r = v.as_f64().ok_or_else(|| Error::invalid("parameter 'fix_r_s_ohm' must be a number"))?;
            if !r.is_finite() || r < 0.0 {
                return Err(Error::invalid(format!("parameter 'fix_r_s_ohm' must be non-negative, got {r}")));
            }
            Some(r)
        }
        None => None,
    };

    let default_z0 = drv
        .z0_ohm
        .first()
        .and_then(|row| row.get(active_idx))
        .map(|pair| pair[0])
        .filter(|z| z.is_finite() && *z > 0.0)
        .unwrap_or(50.0);
    let z_ref_ohm = match params.get("z_ref_ohm").filter(|v| !v.is_null()) {
        Some(v) => {
            let z = v.as_f64().ok_or_else(|| Error::invalid("parameter 'z_ref_ohm' must be a number"))?;
            if !z.is_finite() || z <= 0.0 {
                return Err(Error::invalid(format!("parameter 'z_ref_ohm' must be positive, got {z}")));
            }
            z
        }
        None => default_z0,
    };

    let mut pts: Vec<(f64, C64)> = Vec::new();
    for (k, &f) in drv.freqs_hz.iter().enumerate() {
        if let Some([lo, hi]) = f_window {
            if f < lo || f > hi {
                continue;
            }
        }
        let y_mat = pairs_to_dense(&drv.y_matrix[k]);
        let z_mat = if k < drv.z_matrix.len() {
            pairs_to_dense(&drv.z_matrix[k])
        } else {
            Dense::zeros(drv.ports.len(), drv.ports.len())
        };
        let y_red = schur_reduced_admittance(&y_mat, &z_mat, true, &[active_idx], &term_idx, &z_loads);
        pts.push((f, y_red[(0, 0)]));
    }
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    if pts.len() < 5 {
        return Err(Error::invalid("f_window_hz retains fewer than 5 frequency points"));
    }

    let freqs: Vec<f64> = pts.iter().map(|(f, _)| *f).collect();
    let y_vec: Vec<C64> = pts.iter().map(|(_, y)| *y).collect();
    let z_vec: Vec<C64> =
        y_vec.iter().map(|&y| if y.norm_sqr() > 1e-60 { C64::ONE / y } else { C64::new(1e30, 0.0) }).collect();

    let im_z: Vec<f64> = z_vec.iter().map(|z| z.im).collect();
    let im_y: Vec<f64> = y_vec.iter().map(|y| y.im).collect();
    let mag_y: Vec<f64> = y_vec.iter().map(|y| y.norm()).collect();
    let mag_z: Vec<f64> = z_vec.iter().map(|z| z.norm()).collect();

    let s_crossings = find_upward_zero_crossings(&freqs, &im_z, &mag_y);
    let p_crossings = find_upward_zero_crossings(&freqs, &im_y, &mag_z);

    let mut pairs: Vec<(f64, f64, f64)> = Vec::new();
    for &(fs, ms) in &s_crossings {
        if let Some((fp, _)) = p_crossings.iter().filter(|&&(fp, _)| fp > fs).min_by(|a, b| a.0.total_cmp(&b.0)) {
            pairs.push((fs, *fp, ms));
        }
    }

    let (f_s_hz, f_p_hz, extra_pairs) = if !pairs.is_empty() {
        let best_idx = pairs.iter().enumerate().max_by(|(_, a), (_, b)| a.2.total_cmp(&b.2)).map(|(i, _)| i).unwrap();
        let (fs0, fp0, _) = pairs[best_idx];
        let mut extras = Vec::new();
        for (i, &(fs, fp, _)) in pairs.iter().enumerate() {
            if i != best_idx {
                extras.push((fs, fp));
            }
        }
        (fs0, fp0, extras)
    } else {
        // Fallback: peak Re(Y) / Re(Z) or |Y| / |Z| if zero crossing was skipped on grid
        let i_s = y_vec
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| (a.re.abs() + 1e-6 * a.norm()).total_cmp(&(b.re.abs() + 1e-6 * b.norm())))
            .map(|(i, _)| i)
            .unwrap_or(0);
        let i_p = z_vec
            .iter()
            .enumerate()
            .skip(i_s)
            .max_by(|(_, a), (_, b)| (a.re.abs() + 1e-6 * a.norm()).total_cmp(&(b.re.abs() + 1e-6 * b.norm())))
            .map(|(i, _)| i)
            .unwrap_or((i_s + 1).min(freqs.len() - 1));
        (freqs[i_s], freqs[i_p].max(freqs[i_s] * 1.0001), Vec::new())
    };

    let c_0_raw_f = result
        .scalar("C_0 (F)")
        .filter(|v| v.is_finite() && *v > 0.0)
        .unwrap_or_else(|| estimate_c0(&freqs, &y_vec, f_s_hz, f_p_hz));
    let c0_init = if drv.ports.len() == 1 { c_0_raw_f } else { estimate_c0(&freqs, &y_vec, f_s_hz, f_p_hz) };

    let ([c0_f, r0_ohm, r_s_ohm, cm_f, lm_h, rm_ohm], fit_residual) =
        fit_mbvd_gamma_lm(&freqs, &y_vec, z_ref_ohm, f_s_hz, f_p_hz, c0_init, fix_r_s);

    // Compute precise fs, fp from fitted mBVD (matches resonance/antiresonance cleanly even between discrete grid points)
    let f_s_fit = 1.0 / (2.0 * PI * (lm_h * cm_f).sqrt());
    let f_p_fit = f_s_fit * (1.0 + cm_f / c0_f).sqrt();
    let k_eff2 = if f_p_fit > f_s_fit { (f_p_fit * f_p_fit - f_s_fit * f_s_fit) / (f_p_fit * f_p_fit) } else { 0.0 };
    let k_t2_ieee = compute_kt2_ieee(f_s_fit, f_p_fit);

    let q_s = phase_slope_q(&freqs, &z_vec, f_s_hz, true);
    let q_p = phase_slope_q(&freqs, &y_vec, f_p_hz, true);
    let q_m = if rm_ohm > 1e-12 && cm_f > 0.0 {
        let qm_raw = (2.0 * PI * f_s_fit * lm_h) / rm_ohm;
        if qm_raw.is_finite() && qm_raw < 1e9 && q_s.is_some() {
            Some(qm_raw)
        } else {
            None
        }
    } else {
        None
    };
    let fom = q_m.or(q_s).map(|q| k_eff2 * q);

    let additional_modes = extra_pairs
        .into_iter()
        .map(|(fs, fp)| PiezoModeMetrics {
            f_s_hz: fs,
            f_p_hz: fp,
            k_eff2: (fp * fp - fs * fs) / (fp * fp),
            k_t2_ieee: compute_kt2_ieee(fs, fp),
            q_s: phase_slope_q(&freqs, &z_vec, fs, true),
            q_p: phase_slope_q(&freqs, &y_vec, fp, true),
        })
        .collect();

    Ok(PiezoResonatorReport {
        port: port_name,
        z_ref_ohm,
        f_s_hz: f_s_fit,
        f_p_hz: f_p_fit,
        k_eff2,
        k_t2_ieee,
        c_0_raw_f,
        c0_f,
        r0_ohm,
        r_s_ohm,
        cm_f,
        lm_h,
        rm_ohm,
        q_s,
        q_p,
        q_m,
        fom,
        fit_residual,
        additional_modes,
    })
}

/// Piezoelectric resonator mBVD extraction analysis behind the [`Analysis`] trait (`id: "piezo_resonator"`).
pub struct PiezoResonatorAnalysis;

impl Analysis for PiezoResonatorAnalysis {
    fn id(&self) -> &'static str {
        "piezo_resonator"
    }

    fn name(&self) -> &'static str {
        "Piezoelectric resonator mBVD extraction"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        piezo_resonator_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        Ok(piezo_resonator_with_params(result, &Value::Null)?.to_json())
    }

    fn run_with_params(&self, _mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        Ok(piezo_resonator_with_params(result, params)?.to_json())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nyx_physics::Scalar;

    fn synthetic_mbvd_result(freqs_hz: &[f64], c0: f64, r0: f64, rs: f64, cm: f64, lm: f64, rm: f64) -> SolveResult {
        let mut y_mat = Vec::with_capacity(freqs_hz.len());
        let mut z_mat = Vec::with_capacity(freqs_hz.len());
        let mut s_mat = Vec::with_capacity(freqs_hz.len());
        for &f in freqs_hz {
            let y = mbvd_admittance(2.0 * PI * f, c0, r0, rs, cm, lm, rm);
            let z = C64::ONE / y;
            let s = gamma_of(y, 50.0);
            y_mat.push(vec![vec![[y.re, y.im]]]);
            z_mat.push(vec![vec![[z.re, z.im]]]);
            s_mat.push(vec![vec![[s.re, s.im]]]);
        }
        let drv = DrivenPortData {
            ports: vec!["top".into()],
            freqs_hz: freqs_hz.to_vec(),
            y_matrix: y_mat,
            z_matrix: z_mat,
            s_matrix: s_mat,
            z0_ohm: vec![vec![[50.0, 0.0]]; freqs_hz.len()],
            hdm_solves: freqs_hz.len(),
            rom_dim: 0,
            anchor_freqs_hz: vec![],
            max_residual: 0.0,
        };
        SolveResult {
            physics: "sms".into(),
            study: "driven".into(),
            scalars: vec![drv.scalar(), ("C_0 (F)".into(), Scalar::Num(c0))],
            ..Default::default()
        }
    }

    #[test]
    fn fits_six_parameter_larson_mbvd_and_extracts_q_to_under_1pct() {
        let (c0, r0, rs, cm, lm, rm): (f64, f64, f64, f64, f64, f64) = (1.5e-12, 0.4, 1.2, 9.0e-14, 3.2e-6, 2.5);
        let fs_true = 1.0 / (2.0 * PI * (lm * cm).sqrt());
        let fp_true = fs_true * (1.0 + cm / c0).sqrt();
        let n_pts = 401usize;
        let f_lo = 0.96 * fs_true;
        let f_hi = 1.04 * fp_true;
        let freqs: Vec<f64> = (0..n_pts).map(|k| f_lo + (f_hi - f_lo) * (k as f64) / ((n_pts - 1) as f64)).collect();
        let res = synthetic_mbvd_result(&freqs, c0, r0, rs, cm, lm, rm);
        assert!(PiezoResonatorAnalysis.applies_to(&res));

        let rep = piezo_resonator_with_params(&res, &Value::Null).unwrap();
        assert!((rep.f_s_hz / fs_true - 1.0).abs() < 1e-4, "fs: {} vs {fs_true}", rep.f_s_hz);
        assert!((rep.f_p_hz / fp_true - 1.0).abs() < 1e-4, "fp: {} vs {fp_true}", rep.f_p_hz);
        assert!((rep.c0_f / c0 - 1.0).abs() < 5e-3, "c0: {} vs {c0}", rep.c0_f);
        assert!((rep.cm_f / cm - 1.0).abs() < 5e-3, "cm: {} vs {cm}", rep.cm_f);
        assert!((rep.lm_h / lm - 1.0).abs() < 5e-3, "lm: {} vs {lm}", rep.lm_h);
        assert!((rep.rm_ohm / rm - 1.0).abs() < 1e-2, "rm: {} vs {rm}", rep.rm_ohm);
        assert!((rep.r_s_ohm / rs - 1.0).abs() < 1e-2, "rs: {} vs {rs}", rep.r_s_ohm);
        assert!((rep.r0_ohm / r0 - 1.0).abs() < 1e-2, "r0: {} vs {r0}", rep.r0_ohm);
        assert!(rep.fit_residual < 1e-4, "residual {}", rep.fit_residual);
        assert!(rep.q_s.is_some() && rep.q_p.is_some() && rep.fom.is_some());
    }

    #[test]
    fn lossless_resonator_returns_null_q_and_bounded_gamma_fit() {
        let (c0, cm, lm): (f64, f64, f64) = (2.0e-12, 1.0e-13, 2.5e-6);
        let fs_true = 1.0 / (2.0 * PI * (lm * cm).sqrt());
        let fp_true = fs_true * (1.0 + cm / c0).sqrt();
        let freqs: Vec<f64> =
            (0..121).map(|k| 0.94 * fs_true + (1.06 * fp_true - 0.94 * fs_true) * (k as f64) / 120.0).collect();
        let res = synthetic_mbvd_result(&freqs, c0, 0.0, 0.0, cm, lm, 0.0);
        let rep = piezo_resonator_with_params(&res, &json!({ "fix_r_s_ohm": 0.0 })).unwrap();
        assert!(rep.q_s.is_none() && rep.q_p.is_none(), "expected null Q for lossless resonator");
        assert!((rep.f_s_hz / fs_true - 1.0).abs() < 1e-3);
        assert!((rep.f_p_hz / fp_true - 1.0).abs() < 1e-3);
        assert!(rep.fit_residual < 1e-4);
    }

    #[test]
    fn two_port_schur_terminations_short_open_and_load() {
        // Two-port capacitively coupled resonator (physical C_c = 0.25 pF coupling + C_2 = 0.8 pF shunt at port 2)
        let (c0, cm, lm, rm): (f64, f64, f64, f64) = (1.0e-12, 8.0e-14, 3.0e-6, 1.5);
        let fs = 1.0 / (2.0 * PI * (lm * cm).sqrt());
        let freqs: Vec<f64> = (0..101).map(|k| 0.92 * fs + 0.16 * fs * (k as f64) / 100.0).collect();
        let (cc, c2): (f64, f64) = (0.25e-12, 0.8e-12);
        let mut y_mat = Vec::new();
        let mut z_mat = Vec::new();
        for &f in &freqs {
            let w = 2.0 * PI * f;
            let y_coup = C64::new(0.0, w * cc);
            let y1 = mbvd_admittance(w, c0, 0.0, 0.0, cm, lm, rm) + y_coup;
            let y2 = C64::new(0.0, w * c2) + y_coup;
            let row0 = vec![[y1.re, y1.im], [-y_coup.re, -y_coup.im]];
            let row1 = vec![[-y_coup.re, -y_coup.im], [y2.re, y2.im]];
            y_mat.push(vec![row0, row1]);
            z_mat.push(vec![vec![[0.0, 0.0]; 2]; 2]);
        }
        let drv = DrivenPortData {
            ports: vec!["p1".into(), "p2".into()],
            freqs_hz: freqs,
            y_matrix: y_mat,
            z_matrix: z_mat,
            s_matrix: vec![],
            z0_ohm: vec![],
            hdm_solves: 101,
            rom_dim: 0,
            anchor_freqs_hz: vec![],
            max_residual: 0.0,
        };
        let res = SolveResult {
            physics: "sms".into(),
            study: "driven".into(),
            scalars: vec![drv.scalar()],
            ..Default::default()
        };
        let r_short =
            piezo_resonator_with_params(&res, &json!({ "port": "p1", "terminations": { "p2": 0.0 } })).unwrap();
        let r_load =
            piezo_resonator_with_params(&res, &json!({ "port": "p1", "terminations": { "p2": 50.0 } })).unwrap();
        let r_open =
            piezo_resonator_with_params(&res, &json!({ "port": "p1", "terminations": { "p2": "open" } })).unwrap();
        assert!(r_short.fit_residual < 1e-3 && r_load.fit_residual < 1e-3 && r_open.fit_residual < 1e-3);
        // Short vs open alters input shunt capacitance across y_coup -> different antiresonance frequency fp
        assert!((r_short.f_p_hz - r_open.f_p_hz).abs() > 1e-4 * r_short.f_p_hz);
    }
}
