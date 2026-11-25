// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Canonical mathematical post-processing routines (`Plan Appendix A` & `§5`):
//! - Mode G2 dual-box `1/R` open-boundary Richardson extrapolation (`§2.2.C`);
//! - Charge-basis transmon diagonalization & exact Two-Pole Kubo excited-state admittance `Y_Kubo^(g/e)(omega)`
//!   with Thomas–Reiche–Kuhn oscillator-strength sum-rule verification (`Appendix A.2`);
//! - Algebraic Taubin circle fit + phase-slope resonator fit (`Appendix A.3`);
//! - SQuADDS compositional synthesis (`find_kappa` & claw-transmon `g`, `Appendix A.4`);
//! - Tutorial 9 low-frequency `Im[Y(omega)] / omega` polynomial capacitance extraction (`§3.3.A`);
//! - Ultra-strong damped filter Lorentzian/rational pole extraction (`§3.4.C`);
//! - `C0.1` Schwarz–Christoffel `(c_0, c_1)` `I_edge` thickness calibration & Roache Richardson GCI (`§5.2`).

#![allow(clippy::needless_range_loop)]

use nyx_engine::analysis::epr::{E_CHARGE, PLANCK_H};
use nyx_engine::facade::physics::{DrivenPortData, SolveResult};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Default)]
pub struct Cpx {
    pub re: f64,
    pub im: f64,
}

impl Cpx {
    #[inline]
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    #[inline]
    pub fn from_arr(a: [f64; 2]) -> Self {
        Self { re: a[0], im: a[1] }
    }

    #[inline]
    pub fn to_arr(self) -> [f64; 2] {
        [self.re, self.im]
    }

    #[inline]
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    #[inline]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    #[inline]
    pub fn conj(self) -> Self {
        Self { re: self.re, im: -self.im }
    }

    #[inline]
    pub fn inv(self) -> Self {
        let d = self.norm_sqr().max(1.0e-300);
        Self { re: self.re / d, im: -self.im / d }
    }
}

impl std::ops::Add for Cpx {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self { re: self.re + rhs.re, im: self.im + rhs.im }
    }
}

impl std::ops::Sub for Cpx {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self { re: self.re - rhs.re, im: self.im - rhs.im }
    }
}

impl std::ops::Mul for Cpx {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self { re: self.re * rhs.re - self.im * rhs.im, im: self.re * rhs.im + self.im * rhs.re }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl std::ops::Div for Cpx {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}

impl std::ops::Mul<f64> for Cpx {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: f64) -> Self {
        Self { re: self.re * rhs, im: self.im * rhs }
    }
}

/// Mode G2 dual-box `1/R` open-boundary extrapolation (`Plan §2.2.C`):
/// `C_inf = (R2 * C(R2) - R1 * C(R1)) / (R2 - R1)` with `R_k = 0.5 * L_chip + P_k`.
pub fn extrapolate_dual_box_1_over_r(c_p1: &[Vec<f64>], c_p2: &[Vec<f64>], r1_m: f64, r2_m: f64) -> Vec<Vec<f64>> {
    let n = c_p1.len();
    let denom = (r2_m - r1_m).max(1.0e-18);
    let mut out = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            out[i][j] = (r2_m * c_p2[i][j] - r1_m * c_p1[i][j]) / denom;
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let sym = 0.5 * (out[i][j] + out[j][i]);
            out[i][j] = sym;
            out[j][i] = sym;
        }
    }
    out
}

/// Converts a symmetric Maxwell capacitance matrix `C^M` into mutual/SPICE capacitance matrix `C^mut`.
pub fn mutual_from_maxwell(maxwell: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = maxwell.len();
    let mut mut_mat = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            s += maxwell[i][j];
            if i != j {
                mut_mat[i][j] = -maxwell[i][j];
            }
        }
        mut_mat[i][i] = s;
    }
    mut_mat
}

/// Symmetric dense eigensolver via `femcore::linalg::driver::dense_eigh`, returning ascending `(evals, evecs)`
/// where `evecs[k]` is the `k`-th normalized eigenvector.
pub fn sym_eigh_small(a_in: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a_in.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut flat = Vec::with_capacity(n * n);
    for r in 0..n {
        for c in 0..n {
            flat.push(0.5 * (a_in[r][c] + a_in[c][r]));
        }
    }
    femcore::linalg::driver::dense_eigh(&flat, n).unwrap_or_else(|_| (vec![0.0; n], vec![vec![0.0; n]; n]))
}

/// Checks symmetric Maxwell matrix invariants (`min_eig`, `max_offdiag_pos_ratio`).
pub fn check_maxwell_invariants(maxwell: &[Vec<f64>]) -> (f64, f64) {
    let n = maxwell.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let (evals, _) = sym_eigh_small(maxwell);
    let min_eig = evals.first().copied().unwrap_or(0.0);
    let mut max_pos_ratio = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let denom = (maxwell[i][i].abs() * maxwell[j][j].abs()).sqrt().max(1.0e-30);
            let ratio = maxwell[i][j] / denom;
            if ratio > max_pos_ratio {
                max_pos_ratio = ratio;
            }
        }
    }
    (min_eig, max_pos_ratio)
}

/// Checks driven S-parameter reciprocity `max |S_ij - S_ji|` and lossless unitarity `max ||S^H S - I||_F`.
pub fn check_s_matrix_invariants(s_matrix: &[Vec<Vec<[f64; 2]>>]) -> (f64, f64) {
    let mut max_recip = 0.0_f64;
    let mut max_unit = 0.0_f64;
    for s_f in s_matrix {
        let p = s_f.len();
        for i in 0..p {
            for j in (i + 1)..p {
                let d = Cpx::from_arr(s_f[i][j]) - Cpx::from_arr(s_f[j][i]);
                max_recip = max_recip.max(d.abs());
            }
        }
        // Frobenius upper bound on spectral norm ||S^H S - I||_2 <= ||S^H S - I||_F
        let mut frob_sq = 0.0_f64;
        for i in 0..p {
            for j in 0..p {
                let mut sum = Cpx::new(0.0, 0.0);
                for k in 0..p {
                    sum = sum + Cpx::from_arr(s_f[k][i]).conj() * Cpx::from_arr(s_f[k][j]);
                }
                if i == j {
                    sum.re -= 1.0;
                }
                frob_sq += sum.norm_sqr();
            }
        }
        max_unit = max_unit.max(frob_sq.sqrt());
    }
    (max_recip, max_unit)
}

/// Exact charge-basis transmon spectrum and oscillator strengths (`Plan Appendix A.2`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TransmonKuboModel {
    pub c_env_f: f64,
    pub c_j_f: f64,
    pub c_sigma_f: f64,
    pub e_c_hz: f64,
    pub e_j_hz: f64,
    pub f01_bare_hz: f64,
    pub f12_bare_hz: f64,
    pub alpha_bare_hz: f64,
    pub n01_abs: f64,
    /// `(omega_0m, w_0m)` for `m >= 1`
    pub poles_g: Vec<(f64, f64)>,
    /// `(omega_1m, w_1m)` for `m >= 2`
    pub poles_e_up: Vec<(f64, f64)>,
    pub omega_01: f64,
    pub w_01: f64,
    /// Thomas–Reiche–Kuhn sum `\sum_{m >= 1} w_{0m} == 1`
    pub trk_sum_g: f64,
    /// Thomas–Reiche–Kuhn sum `\sum_{m >= 2} w_{1m} - w_{01} == 1`
    pub trk_sum_e: f64,
}

impl TransmonKuboModel {
    /// Diagonalizes `H / h = 4 E_C (n - n_g)^2 - E_J cos(phi)` in `2*n_charge + 1 = 41` Cooper-pair charge states.
    /// Because `[phi, n] = i` gives `[n, [H/h, n]] = E_J cos(phi)`? No: `[[n, H/h], n]` via kinetic term `4 E_C n^2`:
    /// `[[phi, 4 E_C n^2], phi] = -8 E_C` => TRK canonical sum `sum_{m!=k} sgn(m-k) * (16 E_C |<k|n|m>|^2 / f_{mk}) == 1`
    /// (`|<0|n|1>|^2 = f_01 / (16 E_C)` in the harmonic limit).
    pub fn build(c_env_f: f64, c_j_f: f64, e_j_hz: f64, n_g: f64) -> Self {
        let c_sigma_f = (c_env_f + c_j_f).max(1.0e-18);
        let e_c_hz = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * c_sigma_f);
        let n_max: isize = 20;
        let dim = (2 * n_max + 1) as usize;
        let mut h = vec![vec![0.0; dim]; dim];
        let charges: Vec<f64> = (-n_max..=n_max).map(|n| (n as f64) - n_g).collect();
        for i in 0..dim {
            h[i][i] = 4.0 * e_c_hz * charges[i] * charges[i];
            if i + 1 < dim {
                h[i][i + 1] = -0.5 * e_j_hz;
                h[i + 1][i] = -0.5 * e_j_hz;
            }
        }
        let (evals, evecs) = sym_eigh_small(&h);
        let n_elem = |k: usize, m: usize| -> f64 {
            let mut s = 0.0;
            for i in 0..dim {
                s += charges[i] * evecs[k][i] * evecs[m][i];
            }
            s
        };

        let f01 = (evals[1] - evals[0]).max(1.0);
        let f12 = (evals[2] - evals[1]).max(1.0);
        let alpha = f12 - f01;
        let n01_abs = n_elem(0, 1).abs();
        let w_01 = (16.0 * e_c_hz * n01_abs * n01_abs) / f01;

        let mut poles_g = Vec::new();
        let mut trk_sum_g = 0.0;
        for m in 1..dim.min(18) {
            let f0m = (evals[m] - evals[0]).max(1.0);
            let n0m = n_elem(0, m);
            let w0m = (16.0 * e_c_hz * n0m * n0m) / f0m;
            trk_sum_g += w0m;
            if w0m > 1.0e-15 {
                poles_g.push((2.0 * PI * f0m, w0m));
            }
        }

        let mut poles_e_up = Vec::new();
        let mut trk_up = 0.0;
        for m in 2..dim.min(19) {
            let f1m = (evals[m] - evals[1]).max(1.0);
            let n1m = n_elem(1, m);
            let w1m = (16.0 * e_c_hz * n1m * n1m) / f1m;
            trk_up += w1m;
            if w1m > 1.0e-15 {
                poles_e_up.push((2.0 * PI * f1m, w1m));
            }
        }
        let trk_sum_e = trk_up - w_01;

        Self {
            c_env_f,
            c_j_f,
            c_sigma_f,
            e_c_hz,
            e_j_hz,
            f01_bare_hz: f01,
            f12_bare_hz: f12,
            alpha_bare_hz: alpha,
            n01_abs,
            poles_g,
            poles_e_up,
            omega_01: 2.0 * PI * f01,
            w_01,
            trk_sum_g,
            trk_sum_e,
        }
    }

    /// Intrinsic junction-port Kubo admittance `Y_Kubo^(g/e)(omega) = [Z_Q^(g/e)(omega)]^-1 - j omega C_env` (`+j omega` convention).
    pub fn y_kubo(&self, freq_hz: f64, excited: bool) -> Cpx {
        let omega = 2.0 * PI * freq_hz;
        let om2 = omega * omega;
        let sum_term = if !excited {
            let mut s = 0.0;
            for &(om_0m, w_0m) in &self.poles_g {
                let d = om_0m * om_0m - om2;
                let d_safe = if d.abs() < 1.0e6 { d.signum() * 1.0e6 } else { d };
                s += w_0m / d_safe;
            }
            s
        } else {
            let mut s = 0.0;
            for &(om_1m, w_1m) in &self.poles_e_up {
                let d = om_1m * om_1m - om2;
                let d_safe = if d.abs() < 1.0e6 { d.signum() * 1.0e6 } else { d };
                s += w_1m / d_safe;
            }
            let d01 = self.omega_01 * self.omega_01 - om2;
            let d01_safe = if d01.abs() < 1.0e6 { d01.signum() * 1.0e6 } else { d01 };
            s -= self.w_01 / d01_safe;
            s
        };
        // Z_Q = j * (omega / C_Sigma) * sum_term  =>  1 / Z_Q = -j * C_Sigma / (omega * sum_term)
        let inv_zq_im = if sum_term.abs() > 1.0e-30 { -self.c_sigma_f / (omega * sum_term) } else { 0.0 };
        Cpx::new(0.0, inv_zq_im - omega * self.c_env_f)
    }
}

/// Synthesizes 2-port feedline `DrivenPortData` when junction port `jj_idx` of a 3-port `DrivenPortData` is terminated
/// by `Y_Kubo^(g/e)(omega)` (`Plan Appendix A.2`).
pub fn synthesize_kubo_feedline_2port(
    driven_3p: &DrivenPortData,
    feed_indices: [usize; 2],
    jj_idx: usize,
    z0_feed_ohm: f64,
    kubo: &TransmonKuboModel,
    excited: bool,
) -> DrivenPortData {
    let y0 = 1.0 / z0_feed_ohm.max(1.0e-6);
    let nf = driven_3p.freqs_hz.len();
    let mut s_mat = Vec::with_capacity(nf);
    let mut y_mat = Vec::with_capacity(nf);
    let mut z_mat = Vec::with_capacity(nf);

    for (k, &f_hz) in driven_3p.freqs_hz.iter().enumerate() {
        let y3 = &driven_3p.y_matrix[k];
        let y_jj_tot = Cpx::from_arr(y3[jj_idx][jj_idx]) + kubo.y_kubo(f_hz, excited);
        let inv_yjj = y_jj_tot.inv();

        let mut y_ff = [[Cpx::new(0.0, 0.0); 2]; 2];
        for a in 0..2 {
            let ia = feed_indices[a];
            let y_a_jj = Cpx::from_arr(y3[ia][jj_idx]);
            for b in 0..2 {
                let ib = feed_indices[b];
                let y_jj_b = Cpx::from_arr(y3[jj_idx][ib]);
                y_ff[a][b] = Cpx::from_arr(y3[ia][ib]) - y_a_jj * inv_yjj * y_jj_b;
            }
        }

        // Convert 2x2 Y_ff to S_ff = (Y0 I - Y_ff) (Y0 I + Y_ff)^-1
        let a00 = Cpx::new(y0, 0.0) + y_ff[0][0];
        let a01 = y_ff[0][1];
        let a10 = y_ff[1][0];
        let a11 = Cpx::new(y0, 0.0) + y_ff[1][1];
        let det_a = a00 * a11 - a01 * a10;
        let inv_det_a = det_a.inv();
        let inv_a00 = a11 * inv_det_a;
        let inv_a01 = Cpx::new(-a01.re, -a01.im) * inv_det_a;
        let inv_a10 = Cpx::new(-a10.re, -a10.im) * inv_det_a;
        let inv_a11 = a00 * inv_det_a;

        let n00 = Cpx::new(y0, 0.0) - y_ff[0][0];
        let n01 = Cpx::new(-y_ff[0][1].re, -y_ff[0][1].im);
        let n10 = Cpx::new(-y_ff[1][0].re, -y_ff[1][0].im);
        let n11 = Cpx::new(y0, 0.0) - y_ff[1][1];

        let s00 = n00 * inv_a00 + n01 * inv_a10;
        let s01 = n00 * inv_a01 + n01 * inv_a11;
        let s10 = n10 * inv_a00 + n11 * inv_a10;
        let s11 = n10 * inv_a01 + n11 * inv_a11;

        // Also invert 2x2 Y_ff -> Z_ff
        let det_y = y_ff[0][0] * y_ff[1][1] - y_ff[0][1] * y_ff[1][0];
        let inv_det_y = if det_y.norm_sqr() > 1.0e-60 { det_y.inv() } else { Cpx::new(0.0, 0.0) };
        let z00 = y_ff[1][1] * inv_det_y;
        let z01 = Cpx::new(-y_ff[0][1].re, -y_ff[0][1].im) * inv_det_y;
        let z10 = Cpx::new(-y_ff[1][0].re, -y_ff[1][0].im) * inv_det_y;
        let z11 = y_ff[0][0] * inv_det_y;

        s_mat.push(vec![vec![s00.to_arr(), s01.to_arr()], vec![s10.to_arr(), s11.to_arr()]]);
        y_mat
            .push(vec![vec![y_ff[0][0].to_arr(), y_ff[0][1].to_arr()], vec![y_ff[1][0].to_arr(), y_ff[1][1].to_arr()]]);
        z_mat.push(vec![vec![z00.to_arr(), z01.to_arr()], vec![z10.to_arr(), z11.to_arr()]]);
    }

    // Also note in sym_eigh_small: increase Jacobi max sweeps to 200 for 41x41 tridiagonal convergence!
    DrivenPortData {
        freqs_hz: driven_3p.freqs_hz.clone(),
        ports: vec![driven_3p.ports[feed_indices[0]].clone(), driven_3p.ports[feed_indices[1]].clone()],
        z0_ohm: vec![vec![[z0_feed_ohm, 0.0], [z0_feed_ohm, 0.0]]; nf],
        s_matrix: s_mat,
        z_matrix: z_mat,
        y_matrix: y_mat,
        hdm_solves: driven_3p.hdm_solves,
        rom_dim: driven_3p.rom_dim,
        anchor_freqs_hz: driven_3p.anchor_freqs_hz.clone(),
        max_residual: driven_3p.max_residual,
    }
}

/// Runs Nyx's native Probst–Taubin `S21CircleFit` on a synthesized `DrivenPortData` (`hanger` or `transmission`).
pub fn fit_resonator_s21(dp: &DrivenPortData, branch: &str) -> Result<(f64, f64, f64), String> {
    let mut sr = SolveResult { physics: "ems".into(), study: "driven".into(), order: 1, dim: 3, ..Default::default() };
    sr.scalars.push(dp.scalar());
    let rep = nyx_engine::analysis::network::network_with_params(
        &sr,
        &json!({
            "s21_ports": [dp.ports[0], dp.ports[1]],
            "resonator_branch": branch,
        }),
    )
    .map_err(|e| e.to_string())?;
    let fit = rep.s21_fit.ok_or_else(|| "network s21_fit returned None on synthesized S21 trace".to_string())?;
    Ok((fit.f_r_hz, fit.kappa_hz, fit.q_c))
}

/// SQuADDS Compositional Formulas (`Plan Appendix A.4`):
/// Half-wave `find_kappa`: given uncoupled `f_rough_hz`, CPW `z0_ohm` (`50 Ohm`), coupler `c_tb_f` (`top_to_bottom` mutual)
/// and `c_tg_f` (`top_to_ground` mutual), and `z_load_ohm` (`z0_ohm` for end-coupled, `z0_ohm / 2` for tee-coupled).
pub fn squadds_find_kappa_half_wave(
    f_rough_hz: f64,
    c_tb_f: f64,
    c_tg_f: f64,
    z0_ohm: f64,
    z_load_ohm: f64,
) -> (f64, f64) {
    let omega_rough = 2.0 * PI * f_rough_hz.max(1.0);
    let c_res = PI / (2.0 * omega_rough * z0_ohm.max(1.0));
    let pull = (1.0 + (c_tb_f + c_tg_f).max(0.0) / c_res.max(1.0e-18)).powf(-0.5);
    let omega_r = omega_rough * pull;
    let f_r_hz = omega_r / (2.0 * PI);
    let kappa_rad_s = (omega_r * omega_r * c_tb_f * c_tb_f * z_load_ohm) / (c_res + c_tb_f + c_tg_f).max(1.0e-18);
    let kappa_hz = kappa_rad_s / (2.0 * PI);
    (f_r_hz, kappa_hz)
}

/// SQuADDS Compositional Qubit+Cavity synthesis (`Plan Appendix A.4` & `A.1`):
/// Given `TransmonCross` capacitances (`c_cross_cross_f` Maxwell self, `c_cross_claw_f` mutual, `c_claw_gnd_f` mutual),
/// junction `(c_j_f, e_j_hz)`, resonator `(f_r_hz, quarter_wave, z0_ohm)`, returns `(f01_hz, alpha_hz, g_hz, chi_nyx_hz)`.
#[allow(clippy::too_many_arguments)]
pub fn squadds_compositional_qubit_cavity(
    c_cross_cross_f: f64,
    c_cross_claw_mut_f: f64,
    c_claw_gnd_mut_f: f64,
    c_j_f: f64,
    e_j_hz: f64,
    f_r_hz: f64,
    quarter_wave: bool,
    z0_ohm: f64,
) -> (f64, f64, f64, f64) {
    let c_env = (c_cross_cross_f - c_cross_claw_mut_f * 0.0).max(1.0e-18);
    let kubo = TransmonKuboModel::build(c_env, c_j_f, e_j_hz, 0.0);
    let f01_hz = kubo.f01_bare_hz;
    let alpha_hz = kubo.alpha_bare_hz;

    let k_mode = if quarter_wave { 0.5 } else { 1.0 };
    let omega_r = 2.0 * PI * f_r_hz.max(1.0);
    let c_r_eff = (k_mode * PI) / (2.0 * omega_r * z0_ohm.max(1.0));
    let c_tot_r = (c_r_eff + c_claw_gnd_mut_f.max(0.0)).max(1.0e-18);
    let v_zpf_r = (0.5 * PLANCK_H * f_r_hz / c_tot_r).sqrt();
    let beta = c_cross_claw_mut_f / kubo.c_sigma_f;
    let g_hz = (2.0 * E_CHARGE * beta * v_zpf_r / PLANCK_H) * kubo.n01_abs;

    // Exact non-RWA two-mode dispersive pull matching Plan Appendix A.1 (`chi_Nyx = 2 * chi_Koch < 0` when f01 < f_r):
    let delta = f01_hz - f_r_hz;
    let sigma = f01_hz + f_r_hz;
    let denom_rwa = delta * (delta + alpha_hz);
    let chi_rwa = if denom_rwa.abs() > 1.0e6 { 2.0 * g_hz * g_hz * alpha_hz / denom_rwa } else { 0.0 };
    let denom_cr = sigma * (sigma + alpha_hz);
    let chi_cr = if denom_cr.abs() > 1.0e6 { 2.0 * g_hz * g_hz * alpha_hz / denom_cr } else { 0.0 };
    let chi_nyx_hz = chi_rwa + chi_cr;
    (f01_hz, alpha_hz, g_hz, chi_nyx_hz)
}

/// Tutorial 9 (`S3.1`, `Plan §3.3.A`): Fits `Im[Y_{ab}(omega)] / omega = C_{ab}^{(0)} + c_2 omega^2 + c_4 omega^4`
/// over `[f_min_fit, f_max_fit]` via 3x3 normal equations in `x = (omega / omega_scale)^2`.
pub fn fit_driven_capacitance_omega_zero(dp: &DrivenPortData, f_min_hz: f64, f_max_hz: f64) -> Vec<Vec<f64>> {
    let p = dp.ports.len();
    let mut c0_mat = vec![vec![0.0; p]; p];
    let omega_scale = 2.0 * PI * 2.5e9;

    let selected: Vec<usize> = dp
        .freqs_hz
        .iter()
        .enumerate()
        .filter_map(|(k, &f)| (f >= f_min_hz && f <= f_max_hz && f > 0.0).then_some(k))
        .collect();
    if selected.is_empty() {
        return c0_mat;
    }

    for ia in 0..p {
        for ib in 0..p {
            let mut ata = [[0.0_f64; 3]; 3];
            let mut atb = [0.0_f64; 3];
            for &k in &selected {
                let omega = 2.0 * PI * dp.freqs_hz[k];
                let x = (omega / omega_scale) * (omega / omega_scale);
                let basis = [1.0, x, x * x];
                let y_im = dp.y_matrix[k][ia][ib][1];
                let val = y_im / omega;
                for r in 0..3 {
                    atb[r] += basis[r] * val;
                    for c in 0..3 {
                        ata[r][c] += basis[r] * basis[c];
                    }
                }
            }
            let coeff0 = solve_3x3(ata, atb).map(|sol| sol[0]).unwrap_or_else(|| {
                let k0 = selected[0];
                dp.y_matrix[k0][ia][ib][1] / (2.0 * PI * dp.freqs_hz[k0])
            });
            c0_mat[ia][ib] = coeff0;
        }
    }
    for i in 0..p {
        for j in (i + 1)..p {
            let sym = 0.5 * (c0_mat[i][j] + c0_mat[j][i]);
            c0_mat[i][j] = sym;
            c0_mat[j][i] = sym;
        }
    }
    c0_mat
}

fn solve_3x3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for col in 0..3 {
        let mut piv = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1.0e-30 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let inv_d = 1.0 / a[col][col];
        for r in (col + 1)..3 {
            let f = a[r][col] * inv_d;
            for c in col..3 {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0; 3];
    for r in (0..3).rev() {
        let mut s = b[r];
        for c in (r + 1)..3 {
            s -= a[r][c] * x[c];
        }
        x[r] = s / a[r][r];
    }
    Some(x)
}

/// Extracts a heavily damped resonance pole `(f_f_hz, kappa_f_hz, q_f)` from the Schur-reduced environment
/// impedance/admittance `Re[Z_env(f)]` or `Re[Y_env(f)]` over `[f_lo, f_hi]` via Lorentzian least-squares fit (`Plan §3.4.C`).
pub fn fit_damped_lorentzian_pole(
    freqs_hz: &[f64],
    real_part: &[f64],
    f_lo: f64,
    f_hi: f64,
) -> Option<(f64, f64, f64)> {
    let pts: Vec<(f64, f64)> = freqs_hz
        .iter()
        .copied()
        .zip(real_part.iter().copied())
        .filter(|&(f, r)| f >= f_lo && f <= f_hi && r.is_finite() && r > 0.0)
        .collect();
    if pts.len() < 5 {
        return None;
    }
    let (mut pk_idx, mut pk_val) = (0usize, -1.0_f64);
    for (i, &(_f, r)) in pts.iter().enumerate() {
        if r > pk_val {
            pk_val = r;
            pk_idx = i;
        }
    }
    if pk_val <= 0.0 {
        return None;
    }
    // Linearized inverse Lorentzian near peak: 1 / R(f) = a0 + a1 * df + a2 * df^2
    // Peak at df* = -a1 / (2 a2), R_max = 1 / (a0 - a1^2 / (4 a2)), half-max at a2 * (df - df*)^2 = 1 / R_max => FWHM = 2 / sqrt(a2 * R_max)
    let f0 = pts[pk_idx].0;
    let f_scale = 1.0e9;
    let mut ata = [[0.0_f64; 3]; 3];
    let mut atb = [0.0_f64; 3];
    for &(f, r) in &pts {
        if r >= 0.25 * pk_val {
            let w = (r / pk_val).powi(2);
            let x = (f - f0) / f_scale;
            let basis = [1.0, x, x * x];
            let y = 1.0 / r;
            for ri in 0..3 {
                atb[ri] += w * basis[ri] * y;
                for ci in 0..3 {
                    ata[ri][ci] += w * basis[ri] * basis[ci];
                }
            }
        }
    }
    if let Some([a0, a1, a2]) = solve_3x3(ata, atb) {
        let inv_rmax = a0 - (a1 * a1) / (4.0 * a2.max(1.0e-30));
        if a2 > 0.0 && inv_rmax > 0.0 {
            let df_star = -a1 / (2.0 * a2);
            let f_res = f0 + df_star * f_scale;
            let kappa = 2.0 * f_scale * (inv_rmax / a2).sqrt();
            if f_res > f_lo && f_res < f_hi && kappa > 0.0 {
                return Some((f_res, kappa, f_res / kappa));
            }
        }
    }
    // Fallback half-maximum bracket interpolation
    let half = 0.5 * pk_val;
    let mut f_left = pts.first().unwrap().0;
    for i in (0..pk_idx).rev() {
        if pts[i].1 <= half {
            let t = (half - pts[i].1) / (pts[i + 1].1 - pts[i].1).max(1.0e-30);
            f_left = pts[i].0 + t * (pts[i + 1].0 - pts[i].0);
            break;
        }
    }
    let mut f_right = pts.last().unwrap().0;
    for i in pk_idx..(pts.len() - 1) {
        if pts[i + 1].1 <= half {
            let t = (pts[i].1 - half) / (pts[i].1 - pts[i + 1].1).max(1.0e-30);
            f_right = pts[i].0 + t * (pts[i + 1].0 - pts[i].0);
            break;
        }
    }
    let kappa = (f_right - f_left).max(1.0e3);
    Some((f0, kappa, f0 / kappa))
}

/// Calibrates `(c_0, c_1)` in `\Delta C(t) = t * [c_0 + c_1 ln(s_eff / t)] * \tilde{I}_{edge}` by weighted least squares
/// minimizing relative capacitance error across `(delta_c_exact, s_eff, t, i_edge, c_exact)` samples (`C0.1`).
pub fn fit_sc_thickness_coeffs(samples: &[(f64, f64, f64, f64, f64)]) -> (f64, f64) {
    let mut a00 = 0.0;
    let mut a01 = 0.0;
    let mut a11 = 0.0;
    let mut b0 = 0.0;
    let mut b1 = 0.0;
    for &(dc_exact, s_eff, t, i_edge, c_exact) in samples {
        if !(s_eff > t && t > 0.0 && i_edge > 0.0 && c_exact > 0.0) {
            continue;
        }
        let x0 = (t * i_edge) / c_exact;
        let x1 = x0 * (s_eff / t).ln();
        let y = dc_exact / c_exact;
        a00 += x0 * x0;
        a01 += x0 * x1;
        a11 += x1 * x1;
        b0 += x0 * y;
        b1 += x1 * y;
    }
    let det = a00 * a11 - a01 * a01;
    if det.abs() > 1.0e-30 {
        ((a11 * b0 - a01 * b1) / det, (a00 * b1 - a01 * b0) / det)
    } else {
        (0.082, 0.029)
    }
}

/// Roache Grid Convergence Index (`Plan §5.2`) from 3 ordered grid levels `(y_coarse, y_med, y_fine)` at ratio `r > 1`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GciSummary {
    pub r: f64,
    pub y_coarse: f64,
    pub y_med: f64,
    pub y_fine: f64,
    pub observed_order_raw: f64,
    pub observed_order_clamped: f64,
    pub extrapolated: f64,
    pub gci_fine_pct: f64,
}

pub fn compute_roache_gci(y_coarse: f64, y_med: f64, y_fine: f64, r: f64, p_nominal: f64) -> GciSummary {
    let e32 = y_coarse - y_med;
    let e21 = y_med - y_fine;
    let ratio = if e21.abs() > 1.0e-30 { (e32 / e21).abs() } else { r.powf(p_nominal) };
    let p_raw = if ratio > 1.0 { ratio.ln() / r.ln() } else { p_nominal };
    let p_eff = p_raw.clamp(0.5, 4.0);
    let rp_m1 = (r.powf(p_eff) - 1.0).max(1.0e-12);
    let y_ext = y_fine + (y_fine - y_med) / rp_m1;
    let fs = 1.25;
    let gci_pct = 100.0 * fs * ((y_fine - y_med).abs() / y_fine.abs().max(1.0e-30)) / rp_m1;
    GciSummary {
        r,
        y_coarse,
        y_med,
        y_fine,
        observed_order_raw: p_raw,
        observed_order_clamped: p_eff,
        extrapolated: y_ext,
        gci_fine_pct: gci_pct,
    }
}
