// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified multi-port microwave network, QEP, Black-Box Quantization (BBQ), `S_21` circle-fit, and
//! Touchstone v1.1 analysis (`NetworkAnalysis`, FR-04).
//!
//! Operates on either a 3D cavity [`CircuitRom`] or a driven frequency-domain [`DrivenPortData`] payload:
//! - **Companion QEP (`CircuitRom`)**: Linearizes `(K_tot + i k_0 C_tot - k_0^2 M_tot) z = 0` into a balanced
//!   `2r x 2r` complex eigenproblem solved via [`femcore::linalg::dense::schur`] to extract damped complex
//!   eigenfrequencies `f_re_Hz + i f_im_Hz`, quality factors `q_c`, and linewidths `kappa_Hz = 2 |f_im|`.
//! - **ROM Frequency-Domain Synthesis**: Synthesizes multi-port `Z(omega)`, `Y(omega)`, and `S(omega)` across
//!   `freqs_hz` directly from `CircuitRom`.
//! - **Schur-Complement Black-Box Quantization (`bbq_modes`)**: Partitions ports into active junction ports `A`
//!   and terminated ports `T`, forms `Y_red(omega) = Y_AA - Y_AT (Y_TT + Y_load)^-1 Y_TA`, locates Foster
//!   parallel-resonance zeros `Im[Y_JJ(omega_m)] = 0`, and extracts `C_{Sigma,m}`, `L_{eff,m}`, `kappa_Hz`,
//!   `q_m`, and `t1_purcell_s`.
//! - **Taubin Algebraic Circle Fit (`s21_fit`, Probst et al. 2015)**: Fits `S_{out,in}(f)` on the complex plane
//!   via centered `3 x 3` Taubin algebraic circle fitting + Gauss–Newton arctan phase fitting for `"hanger"`,
//!   `"transmission"`, and `"reflection"` resonator topologies (`f_r_Hz`, `q_l`, `q_c`, `q_i`, `kappa_Hz`, `phi_0_rad`).
//! - **Touchstone v1.1 (`.sNp`)**: Formats multi-port scattering parameters in `"RI"`, `"MA"`, or `"DB"` format.
#![allow(clippy::needless_range_loop)]

use crate::epr::{assemble_rom_damping, is_nonlinear_port, jj_inv_l_tot, resolve_ports_with_overrides, C0, EPS0, MU0};
use crate::Analysis;
use femcore::linalg::dense::{schur, solve_small, Dense};
use femcore::linalg::driver::{dense_eigh, dense_gen_eigh};
use femcore::linalg::scalar::{Scalar, C64};
use nyx_physics::{CircuitRom, DrivenPortData, Error, Mesh, RomPort, SolveResult};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::fmt::Write as _;
use std::time::Instant;

/// Complex eigenmode extracted from the `2r x 2r` damped `CircuitRom` quadratic eigenvalue problem.
#[derive(Clone, Debug, PartialEq)]
pub struct QepMode {
    pub f_re_hz: f64,
    pub f_im_hz: f64,
    pub q_c: f64,
    pub kappa_hz: f64,
}

impl QepMode {
    pub fn to_json(&self) -> Value {
        json!({
            "f_re_Hz": self.f_re_hz,
            "f_im_Hz": self.f_im_hz,
            "q_c": self.q_c,
            "kappa_Hz": self.kappa_hz,
        })
    }
}

/// Foster zero-crossing mode extracted by Schur-reduced Black-Box Quantization (Nigg et al. 2012).
#[derive(Clone, Debug, PartialEq)]
pub struct BbqMode {
    pub port: String,
    pub f_m_hz: f64,
    pub c_sigma_f: f64,
    pub l_eff_h: f64,
    pub kappa_hz: f64,
    pub q_m: f64,
    pub t1_purcell_s: f64,
}

impl BbqMode {
    pub fn to_json(&self) -> Value {
        json!({
            "port": self.port,
            "f_m_Hz": self.f_m_hz,
            "C_Sigma_F": self.c_sigma_f,
            "L_eff_H": self.l_eff_h,
            "kappa_Hz": self.kappa_hz,
            "q_m": self.q_m,
            "t1_purcell_s": self.t1_purcell_s,
        })
    }
}

/// Algebraic Taubin circle + phase-vs-frequency resonator fit (Probst et al. 2015).
#[derive(Clone, Debug, PartialEq)]
pub struct S21CircleFit {
    pub ports: [String; 2],
    pub resonator_branch: String,
    pub f_r_hz: f64,
    pub q_l: f64,
    pub q_c: f64,
    pub q_i: f64,
    pub kappa_hz: f64,
    pub phi_0_rad: f64,
    pub circle_center: [f64; 2],
    pub circle_radius: f64,
}

impl S21CircleFit {
    pub fn to_json(&self) -> Value {
        json!({
            "ports": self.ports,
            "resonator_branch": self.resonator_branch,
            "f_r_Hz": self.f_r_hz,
            "q_l": self.q_l,
            "q_c": self.q_c,
            "q_i": self.q_i,
            "kappa_Hz": self.kappa_hz,
            "phi_0_rad": self.phi_0_rad,
            "circle_center": self.circle_center,
            "circle_radius": self.circle_radius,
        })
    }
}

/// Typed output of [`network`] and [`network_with_params`].
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkReport {
    pub ports: Vec<String>,
    pub freqs_hz: Vec<f64>,
    pub z0_ohm: Vec<Vec<[f64; 2]>>,
    pub s_matrix: Vec<Vec<Vec<[f64; 2]>>>,
    pub z_matrix: Vec<Vec<Vec<[f64; 2]>>>,
    pub y_matrix: Vec<Vec<Vec<[f64; 2]>>>,
    pub qep_modes: Vec<QepMode>,
    pub bbq_modes: Vec<BbqMode>,
    pub s21_fit: Option<S21CircleFit>,
    pub touchstone_snp: Option<String>,
    pub elapsed_s: f64,
}

impl NetworkReport {
    pub fn to_json(&self) -> Value {
        let mut obj = json!({
            "ports": self.ports,
            "freqs_Hz": self.freqs_hz,
            "z0_ohm": self.z0_ohm,
            "s_matrix": self.s_matrix,
            "z_matrix": self.z_matrix,
            "y_matrix": self.y_matrix,
            "qep_modes": self.qep_modes.iter().map(QepMode::to_json).collect::<Vec<_>>(),
            "bbq_modes": self.bbq_modes.iter().map(BbqMode::to_json).collect::<Vec<_>>(),
            "elapsed_s": self.elapsed_s,
        });
        if let Some(ref fit) = self.s21_fit {
            obj["s21_fit"] = fit.to_json();
        }
        if let Some(ref snp) = self.touchstone_snp {
            obj["touchstone_snp"] = json!(snp);
        }
        obj
    }
}

/// True when `result` carries a [`CircuitRom`] or [`DrivenPortData`] payload.
pub fn network_applies(result: &SolveResult) -> bool {
    CircuitRom::of(result).is_some() || DrivenPortData::of(result).is_some()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ResonatorBranch {
    Hanger,
    Transmission,
    Reflection,
}

impl ResonatorBranch {
    fn parse(s: &str) -> Result<Self, Error> {
        match s {
            "hanger" => Ok(Self::Hanger),
            "transmission" => Ok(Self::Transmission),
            "reflection" => Ok(Self::Reflection),
            other => Err(Error::invalid(format!(
                "resonator_branch must be 'hanger', 'transmission', or 'reflection', not '{other}'"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Hanger => "hanger",
            Self::Transmission => "transmission",
            Self::Reflection => "reflection",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TouchstoneFormat {
    Ri,
    Ma,
    Db,
}

impl TouchstoneFormat {
    fn parse(s: &str) -> Result<Self, Error> {
        match s.to_ascii_uppercase().as_str() {
            "RI" => Ok(Self::Ri),
            "MA" => Ok(Self::Ma),
            "DB" => Ok(Self::Db),
            _ => Err(Error::invalid(format!("touchstone_format must be 'RI', 'MA', or 'DB', not '{s}'"))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Ri => "RI",
            Self::Ma => "MA",
            Self::Db => "DB",
        }
    }
}

struct NetworkParams {
    terminations_ohm: BTreeMap<String, f64>,
    active_ports: Option<Vec<String>>,
    s21_ports: Option<[String; 2]>,
    resonator_branch: ResonatorBranch,
    has_explicit_branch: bool,
    freqs_hz: Option<Vec<f64>>,
    touchstone_format: Option<TouchstoneFormat>,
}

impl NetworkParams {
    fn from_json(params: &Value) -> Result<Self, Error> {
        if params.is_null() {
            return Ok(Self {
                terminations_ohm: BTreeMap::new(),
                active_ports: None,
                s21_ports: None,
                resonator_branch: ResonatorBranch::Hanger,
                has_explicit_branch: false,
                freqs_hz: None,
                touchstone_format: None,
            });
        }
        let map = params
            .as_object()
            .ok_or_else(|| Error::invalid("network analysis params must be a JSON object or null"))?;

        let mut terminations_ohm = BTreeMap::new();
        if let Some(v) = map.get("terminations_ohm").filter(|v| !v.is_null()) {
            let obj =
                v.as_object().ok_or_else(|| Error::invalid("terminations_ohm must be a port-keyed JSON object"))?;
            for (k, val) in obj {
                let r = val.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                    Error::invalid(format!("terminations_ohm['{k}'] must be a positive finite impedance in Ohm"))
                })?;
                terminations_ohm.insert(k.clone(), r);
            }
        }

        let active_ports = match map.get("active_ports") {
            None | Some(Value::Null) => None,
            Some(Value::Array(arr)) => {
                let mut out = Vec::with_capacity(arr.len());
                for (i, item) in arr.iter().enumerate() {
                    let s = item
                        .as_str()
                        .ok_or_else(|| Error::invalid(format!("active_ports[{i}] must be a port name string")))?;
                    out.push(s.to_string());
                }
                Some(out)
            }
            _ => return Err(Error::invalid("active_ports must be an array of port name strings")),
        };

        let s21_ports = match map.get("s21_ports") {
            None | Some(Value::Null) => None,
            Some(Value::Array(arr)) => {
                if arr.len() != 2 {
                    return Err(Error::invalid("s21_ports must be a 2-element array [input_port, output_port]"));
                }
                let p0 = arr[0]
                    .as_str()
                    .ok_or_else(|| Error::invalid("s21_ports[0] must be a port name string"))?
                    .to_string();
                let p1 = arr[1]
                    .as_str()
                    .ok_or_else(|| Error::invalid("s21_ports[1] must be a port name string"))?
                    .to_string();
                Some([p0, p1])
            }
            _ => return Err(Error::invalid("s21_ports must be a 2-element array [input_port, output_port]")),
        };

        let (resonator_branch, has_explicit_branch) = match map.get("resonator_branch") {
            None | Some(Value::Null) => (ResonatorBranch::Hanger, false),
            Some(Value::String(s)) => (ResonatorBranch::parse(s)?, true),
            _ => return Err(Error::invalid("resonator_branch must be 'hanger', 'transmission', or 'reflection'")),
        };

        let freqs_hz = match map.get("freqs_hz") {
            None | Some(Value::Null) => None,
            Some(Value::Array(arr)) => {
                if arr.is_empty() {
                    return Err(Error::invalid("freqs_hz must not be empty"));
                }
                let mut f_vec = Vec::with_capacity(arr.len());
                for (i, v) in arr.iter().enumerate() {
                    let f = v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                        Error::invalid(format!("freqs_hz[{i}] must be a positive finite frequency in Hz"))
                    })?;
                    if let Some(&prev) = f_vec.last() {
                        if f <= prev {
                            return Err(Error::invalid("freqs_hz must be strictly ascending"));
                        }
                    }
                    f_vec.push(f);
                }
                Some(f_vec)
            }
            _ => return Err(Error::invalid("freqs_hz must be an array of positive numbers")),
        };

        let touchstone_format = match map.get("touchstone_format") {
            None | Some(Value::Null) => None,
            Some(Value::String(s)) => Some(TouchstoneFormat::parse(s)?),
            _ => return Err(Error::invalid("touchstone_format must be 'RI', 'MA', or 'DB'")),
        };

        Ok(Self {
            terminations_ohm,
            active_ports,
            s21_ports,
            resonator_branch,
            has_explicit_branch,
            freqs_hz,
            touchstone_format,
        })
    }
}

// ------------------------------------------------------------------ 1. Companion QEP & CircuitRom S/Z/Y Synthesis

fn assemble_rom_km(rom: &CircuitRom, ports: &[RomPort]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let r = rom.k_r.len();
    let mut k_tot = rom.k_r.clone();
    let mut m_tot = rom.m_r.clone();
    for (p, port) in ports.iter().enumerate() {
        if !is_nonlinear_port(port) {
            continue;
        }
        let inv_l = jj_inv_l_tot(port.l_j_h.unwrap_or(0.0), port.e_l_hz.unwrap_or(0.0));
        let ak = MU0 * inv_l;
        let am = port.c_j_f / EPS0;
        for i in 0..r {
            let bi = rom.b_r[i].get(p).copied().unwrap_or(0.0);
            if bi == 0.0 {
                continue;
            }
            for j in 0..r {
                let bj = rom.b_r[j].get(p).copied().unwrap_or(0.0);
                k_tot[i][j] += ak * bi * bj;
                m_tot[i][j] += am * bi * bj;
            }
        }
    }
    (k_tot, m_tot)
}

/// Solves the `2r x 2r` companion QEP `(K_tot + i k_0 C_tot - k_0^2 M_tot) z = 0` via `femcore::linalg::dense::schur`.
fn solve_rom_qep(rom: &CircuitRom, ports: &[RomPort], c_tot: &[Vec<f64>]) -> Result<Vec<QepMode>, Error> {
    let r = rom.k_r.len();
    if r == 0 {
        return Ok(Vec::new());
    }
    let (k_tot, m_tot) = assemble_rom_km(rom, ports);

    let mut m_flat = Vec::with_capacity(r * r);
    for row in &m_tot {
        m_flat.extend_from_slice(row);
    }
    let (m_eigs, m_vecs) = dense_eigh(&m_flat, r).map_err(Error::from_kernel)?;
    let m_max = m_eigs.iter().copied().fold(0.0_f64, f64::max);
    if m_max <= 0.0 {
        return Err(Error::invalid("CircuitRom mass matrix M_tot is not positive definite"));
    }

    // Build M_tot^{-1/2} = V diag(lam^{-1/2}) V^T
    let mut m_inv_sqrt = vec![vec![0.0; r]; r];
    for (k, &lam) in m_eigs.iter().enumerate() {
        if lam <= 1.0e-14 * m_max {
            return Err(Error::invalid("CircuitRom mass matrix M_tot is numerically singular"));
        }
        let s = 1.0 / lam.sqrt();
        let vk = &m_vecs[k];
        for i in 0..r {
            for j in 0..r {
                m_inv_sqrt[i][j] += s * vk[i] * vk[j];
            }
        }
    }

    let cong = |a: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let mut am = vec![vec![0.0; r]; r];
        for i in 0..r {
            for j in 0..r {
                let mut s = 0.0;
                for k in 0..r {
                    s += a[i][k] * m_inv_sqrt[k][j];
                }
                am[i][j] = s;
            }
        }
        let mut out = vec![vec![0.0; r]; r];
        for i in 0..r {
            for j in i..r {
                let mut s = 0.0;
                for k in 0..r {
                    s += m_inv_sqrt[k][i] * am[k][j];
                }
                out[i][j] = s;
                out[j][i] = s;
            }
        }
        out
    };

    let k_tilde = cong(&k_tot);
    let c_tilde = if c_tot.len() == r { cong(c_tot) } else { vec![vec![0.0; r]; r] };

    let mut k_frob2 = 0.0;
    for row in &k_tilde {
        for &v in row {
            k_frob2 += v * v;
        }
    }
    let s0 = (k_frob2.sqrt() / (r as f64).sqrt()).sqrt().max(1.0);
    let inv_s0 = 1.0 / s0;
    let inv_s0_sq = inv_s0 * inv_s0;

    let mut comp = Dense::<C64>::zeros(2 * r, 2 * r);
    for i in 0..r {
        comp[(i, r + i)] = C64::ONE;
        for j in 0..r {
            comp[(r + i, j)] = C64::new(k_tilde[i][j] * inv_s0_sq, 0.0);
            comp[(r + i, r + j)] = C64::new(0.0, c_tilde[i][j] * inv_s0);
        }
    }

    let (t_mat, _) = schur(&comp).map_err(|e| Error::from_kernel(e.to_string()))?;
    let k0_to_hz = C0 / (2.0 * PI);
    let mut cands: Vec<QepMode> = Vec::new();
    for i in 0..(2 * r) {
        let k0 = t_mat[(i, i)].scale(s0);
        let f_re = k0.re * k0_to_hz;
        let f_im = k0.im * k0_to_hz;
        if f_re > 1.0e3 {
            let kappa_hz = 2.0 * f_im.abs();
            let q_c = if kappa_hz > 1.0e-15 * f_re { f_re / kappa_hz } else { f64::INFINITY };
            cands.push(QepMode { f_re_hz: f_re, f_im_hz: f_im, q_c, kappa_hz });
        }
    }
    cands.sort_by(|a, b| a.f_re_hz.partial_cmp(&b.f_re_hz).unwrap());

    if !rom.modes.is_empty() && cands.len() > rom.modes.len() {
        let mut used = vec![false; cands.len()];
        let mut picked = Vec::with_capacity(rom.modes.len());
        for row in &rom.modes {
            let mut best = None;
            let mut best_d = f64::INFINITY;
            for (k, cm) in cands.iter().enumerate() {
                if !used[k] {
                    let d = (cm.f_re_hz - row.freq_hz).abs() / row.freq_hz.max(1.0);
                    if d < best_d {
                        best_d = d;
                        best = Some(k);
                    }
                }
            }
            if let Some(k) = best {
                used[k] = true;
                picked.push(cands[k].clone());
            }
        }
        picked.sort_by(|a, b| a.f_re_hz.partial_cmp(&b.f_re_hz).unwrap());
        Ok(picked)
    } else if rom.r_static > 0 && rom.modes.is_empty() {
        Ok(cands.into_iter().filter(|m| m.f_re_hz < 50.0e9).collect())
    } else {
        Ok(cands)
    }
}

/// Inverts a `p x p` complex matrix `a` with tiny diagonal regularization fallback if singular.
fn invert_complex_matrix(a: &Dense<C64>) -> Dense<C64> {
    let p = a.rows();
    if p == 0 {
        return Dense::zeros(0, 0);
    }
    let scale = a.frobenius().max(1.0e-30);
    let mut reg_a = a.clone();
    let mut out = Dense::<C64>::zeros(p, p);
    for attempt in 0..2 {
        let mut ok = true;
        for col in 0..p {
            let mut rhs = vec![C64::ZERO; p];
            rhs[col] = C64::ONE;
            if let Some(sol) = solve_small(&reg_a, &rhs) {
                for r in 0..p {
                    out[(r, col)] = sol[r];
                }
            } else {
                ok = false;
                break;
            }
        }
        if ok {
            return out;
        }
        if attempt == 0 {
            for d in 0..p {
                reg_a[(d, d)] += C64::new(1.0e-14 * scale, 0.0);
            }
        }
    }
    out
}

/// Converts a `p x p` impedance matrix `Z` and positive reference impedances `z0_ports` into `S = (Z_n - I)(Z_n + I)^-1`.
fn z_to_s_matrix(z_mat: &Dense<C64>, z0_ports: &[f64]) -> Dense<C64> {
    let p = z_mat.rows();
    if p == 0 {
        return Dense::zeros(0, 0);
    }
    let inv_sqrt_z0: Vec<f64> = z0_ports.iter().map(|&z| 1.0 / z.max(1.0e-30).sqrt()).collect();
    let mut zn_minus = Dense::<C64>::zeros(p, p);
    let mut zn_plus = Dense::<C64>::zeros(p, p);
    for i in 0..p {
        for j in 0..p {
            let zn_ij = z_mat[(i, j)].scale(inv_sqrt_z0[i] * inv_sqrt_z0[j]);
            let d = if i == j { C64::ONE } else { C64::ZERO };
            zn_minus[(i, j)] = zn_ij - d;
            zn_plus[(i, j)] = zn_ij + d;
        }
    }
    let zn_plus_inv = invert_complex_matrix(&zn_plus);
    zn_minus.matmul(&zn_plus_inv)
}

/// Evaluates the open-circuit impedance `Z(omega)` (`P x P`) and short-circuit admittance `Y(omega)` (`P x P`)
/// of `rom` at angular frequency `omega` (rad/s) using passive internal loss `c_internal` (excluding external port loads).
fn eval_rom_zy_at_omega(
    rom: &CircuitRom,
    k_tot: &[Vec<f64>],
    m_tot: &[Vec<f64>],
    c_internal: &[Vec<f64>],
    omega: f64,
) -> Option<(Dense<C64>, Dense<C64>)> {
    let r = rom.k_r.len();
    let p_ports = rom.ports.len();
    if r == 0 || p_ports == 0 {
        return None;
    }
    let k0 = omega / C0;
    let mut a_mat = Dense::<C64>::zeros(r, r);
    for i in 0..r {
        for j in 0..r {
            let re = k_tot[i][j] - k0 * k0 * m_tot[i][j];
            let im = if c_internal.len() == r { k0 * c_internal[i][j] } else { 0.0 };
            a_mat[(i, j)] = C64::new(re, im);
        }
    }

    let mut z_mat = Dense::<C64>::zeros(p_ports, p_ports);
    let mut x_cols: Vec<Vec<C64>> = Vec::with_capacity(p_ports);
    for q in 0..p_ports {
        let rhs: Vec<C64> = (0..r).map(|i| C64::new(rom.b_r[i].get(q).copied().unwrap_or(0.0), 0.0)).collect();
        let sol = if let Some(s) = solve_small(&a_mat, &rhs) {
            s
        } else {
            let mut reg = a_mat.clone();
            let sc = reg.frobenius().max(1.0e-30);
            for d in 0..r {
                reg[(d, d)] += C64::new(0.0, 1.0e-14 * sc);
            }
            solve_small(&reg, &rhs)?
        };
        x_cols.push(sol);
    }

    let pref = C64::new(0.0, omega * MU0);
    for p in 0..p_ports {
        for q in 0..p_ports {
            let mut acc = C64::ZERO;
            for i in 0..r {
                let bp = rom.b_r[i].get(p).copied().unwrap_or(0.0);
                if bp != 0.0 {
                    acc += x_cols[q][i].scale(bp);
                }
            }
            z_mat[(p, q)] = pref * acc;
        }
    }
    let y_mat = invert_complex_matrix(&z_mat);
    Some((z_mat, y_mat))
}

fn dense_to_pairs(m: &Dense<C64>) -> Vec<Vec<[f64; 2]>> {
    (0..m.rows()).map(|i| (0..m.cols()).map(|j| [m[(i, j)].re, m[(i, j)].im]).collect()).collect()
}

fn pairs_to_dense(m: &[Vec<[f64; 2]>]) -> Dense<C64> {
    let r = m.len();
    let c = m.first().map_or(0, Vec::len);
    Dense::from_fn(r, c, |i, j| C64::new(m[i][j][0], m[i][j][1]))
}

// ------------------------------------------------------------------ 2. Schur-Complement Black-Box Quantization (BBQ)

/// Computes `Y_red = Y_AA - Y_AT (Y_TT + Y_load)^-1 Y_TA`, using the algebraically identical dual
/// impedance Schur form `(Z_AA - Z_AT (Z_TT + Z_load)^-1 Z_TA)^-1` whenever `Z` is rank-deficient
/// (`P > r`, where `Y = Z^-1` is regularized) or near an admittance pole.
fn schur_reduced_admittance(
    y_full: &Dense<C64>,
    z_full: &Dense<C64>,
    y_is_exact: bool,
    active_idx: &[usize],
    term_idx: &[usize],
    z_load_ohm: &[f64],
) -> Dense<C64> {
    let na = active_idx.len();
    let nt = term_idx.len();
    if na == 0 {
        return Dense::zeros(0, 0);
    }
    if nt == 0 {
        if y_is_exact {
            return Dense::from_fn(na, na, |i, j| y_full[(active_idx[i], active_idx[j])]);
        }
        let z_aa = Dense::from_fn(na, na, |i, j| z_full[(active_idx[i], active_idx[j])]);
        return invert_complex_matrix(&z_aa);
    }

    if y_is_exact && y_full.is_finite() {
        let y_aa = Dense::from_fn(na, na, |i, j| y_full[(active_idx[i], active_idx[j])]);
        let y_at = Dense::from_fn(na, nt, |i, j| y_full[(active_idx[i], term_idx[j])]);
        let y_ta = Dense::from_fn(nt, na, |i, j| y_full[(term_idx[i], active_idx[j])]);
        let mut y_tt = Dense::from_fn(nt, nt, |i, j| y_full[(term_idx[i], term_idx[j])]);
        for (t, &zl) in z_load_ohm.iter().enumerate() {
            if zl.is_finite() && zl > 0.0 {
                y_tt[(t, t)] += C64::new(1.0 / zl, 0.0);
            }
        }
        let mut x_ta = Dense::<C64>::zeros(nt, na);
        let mut y_ok = true;
        for col in 0..na {
            let rhs: Vec<C64> = (0..nt).map(|r| y_ta[(r, col)]).collect();
            if let Some(sol) = solve_small(&y_tt, &rhs) {
                for r in 0..nt {
                    x_ta[(r, col)] = sol[r];
                }
            } else {
                y_ok = false;
                break;
            }
        }
        if y_ok {
            let corr = y_at.matmul(&x_ta);
            return Dense::from_fn(na, na, |i, j| y_aa[(i, j)] - corr[(i, j)]);
        }
    }

    // Exact Frobenius-Schur / Woodbury dual form on Z (immune to P > r rank deficiency of Z)
    let z_aa = Dense::from_fn(na, na, |i, j| z_full[(active_idx[i], active_idx[j])]);
    let z_at = Dense::from_fn(na, nt, |i, j| z_full[(active_idx[i], term_idx[j])]);
    let z_ta = Dense::from_fn(nt, na, |i, j| z_full[(term_idx[i], active_idx[j])]);
    let mut z_tt = Dense::from_fn(nt, nt, |i, j| z_full[(term_idx[i], term_idx[j])]);
    for (t, &zl) in z_load_ohm.iter().enumerate() {
        z_tt[(t, t)] += C64::new(zl, 0.0);
    }
    let z_tt_inv = invert_complex_matrix(&z_tt);
    let z_corr = z_at.matmul(&z_tt_inv.matmul(&z_ta));
    let z_red = Dense::from_fn(na, na, |i, j| z_aa[(i, j)] - z_corr[(i, j)]);
    invert_complex_matrix(&z_red)
}

/// Evaluates a 4-point Lagrange cubic polynomial and its derivative at `x`.
fn eval_lagrange4_and_deriv(xs: [f64; 4], ys: [f64; 4], x: f64) -> (f64, f64) {
    let mut val = 0.0;
    let mut der = 0.0;
    for i in 0..4 {
        let mut li = 1.0;
        let mut denom = 1.0;
        for j in 0..4 {
            if i != j {
                li *= x - xs[j];
                denom *= xs[i] - xs[j];
            }
        }
        let w = ys[i] / denom;
        val += w * li;
        let mut sum_prod = 0.0;
        for k in 0..4 {
            if k == i {
                continue;
            }
            let mut p = 1.0;
            for j in 0..4 {
                if j != i && j != k {
                    p *= x - xs[j];
                }
            }
            sum_prod += p;
        }
        der += w * sum_prod;
    }
    (val, der)
}

#[allow(clippy::too_many_arguments)]
fn extract_bbq_modes(
    ports: &[String],
    freqs_hz: &[f64],
    y_mats: &[Dense<C64>],
    z_mats: &[Dense<C64>],
    z0_grid: &[Vec<f64>],
    active_idx: &[usize],
    term_idx: &[usize],
    terminations_ohm: &BTreeMap<String, f64>,
    qep_hints_hz: &[f64],
    exact_rom_eval: Option<&dyn Fn(f64) -> Option<Dense<C64>>>,
) -> Vec<BbqMode> {
    let n_pts = freqs_hz.len();
    if n_pts < 2 || active_idx.is_empty() {
        return Vec::new();
    }

    let y_is_exact = exact_rom_eval.is_none();
    let y_red_grid: Vec<Dense<C64>> = (0..n_pts)
        .map(|k| {
            let z_loads: Vec<f64> = term_idx
                .iter()
                .map(|&ti| {
                    terminations_ohm
                        .get(&ports[ti])
                        .copied()
                        .unwrap_or_else(|| z0_grid[k].get(ti).copied().unwrap_or(50.0))
                })
                .collect();
            schur_reduced_admittance(&y_mats[k], &z_mats[k], y_is_exact, active_idx, term_idx, &z_loads)
        })
        .collect();

    let refine_rom_bracket = |eval_fn: &dyn Fn(f64) -> Option<Dense<C64>>,
                              a_slot: usize,
                              port_name: &str,
                              mut lo: f64,
                              mut hi: f64,
                              mut flo: f64,
                              mut fhi: f64|
     -> Option<BbqMode> {
        let mut wm = lo - flo * (hi - lo) / (fhi - flo).max(1e-300);
        for _ in 0..20 {
            wm = wm.clamp(lo, hi);
            let ym_mat = eval_fn(wm)?;
            let fm = ym_mat[(a_slot, a_slot)].im;
            if fm.abs() <= 1.0e-15 * (flo.abs() + fhi.abs()).max(1e-30) || (hi - lo) <= 1.0e-13 * wm {
                break;
            }
            if fm > 0.0 {
                hi = wm;
                fhi = fm;
            } else {
                lo = wm;
                flo = fm;
            }
            let sec = lo - flo * (hi - lo) / (fhi - flo);
            wm = if sec > lo && sec < hi { sec } else { 0.5 * (lo + hi) };
        }
        let dw = (1.0e-7 * wm).max(1.0);
        let (ym, yp, yn) = (eval_fn(wm)?, eval_fn(wm + dw)?, eval_fn(wm - dw)?);
        let y_center = ym[(a_slot, a_slot)];
        let d_im_dw = (yp[(a_slot, a_slot)].im - yn[(a_slot, a_slot)].im) / (2.0 * dw);
        if d_im_dw <= 0.0 {
            return None;
        }
        let c_sigma_f = 0.5 * d_im_dw;
        let f_m_hz = wm / (2.0 * PI);
        let l_eff_h = 1.0 / (wm * wm * c_sigma_f);
        let re_y = y_center.re.max(0.0);
        let kappa_rad_s = re_y / c_sigma_f;
        let kappa_hz = kappa_rad_s / (2.0 * PI);
        let q_m = if kappa_hz > 1.0e-15 * f_m_hz { f_m_hz / kappa_hz } else { f64::INFINITY };
        let t1_purcell_s = if kappa_rad_s > 1.0e-30 { 1.0 / kappa_rad_s } else { f64::INFINITY };
        Some(BbqMode { port: port_name.to_string(), f_m_hz, c_sigma_f, l_eff_h, kappa_hz, q_m, t1_purcell_s })
    };

    let mut modes: Vec<BbqMode> = Vec::new();
    let push_dedup = |modes: &mut Vec<BbqMode>, cand: BbqMode| {
        if !modes.iter().any(|m| m.port == cand.port && (m.f_m_hz - cand.f_m_hz).abs() <= 1.0e-5 * cand.f_m_hz) {
            modes.push(cand);
        }
    };

    for (a_slot, &orig_p) in active_idx.iter().enumerate() {
        let port_name = &ports[orig_p];
        for k in 0..(n_pts - 1) {
            let y0 = y_red_grid[k][(a_slot, a_slot)];
            let y1 = y_red_grid[k + 1][(a_slot, a_slot)];
            if !(y0.is_finite() && y1.is_finite()) {
                continue;
            }
            if !(y0.im <= 0.0 && y1.im > 0.0) {
                continue;
            }
            let w0 = 2.0 * PI * freqs_hz[k];
            let w1 = 2.0 * PI * freqs_hz[k + 1];
            if w1 <= w0 {
                continue;
            }

            if let Some(eval_fn) = exact_rom_eval {
                if let Some(m) = refine_rom_bracket(eval_fn, a_slot, port_name, w0, w1, y0.im, y1.im) {
                    push_dedup(&mut modes, m);
                }
            } else {
                // Discrete DrivenPortData interpolation (4-point cubic when n_pts >= 4, linear otherwise)
                let mut wm = w0 - y0.im * (w1 - w0) / (y1.im - y0.im).max(1e-300);
                let (re_y, d_im_dw) = if n_pts >= 4 {
                    let i0 = k.saturating_sub(1).min(n_pts - 4);
                    let ws = [
                        2.0 * PI * freqs_hz[i0],
                        2.0 * PI * freqs_hz[i0 + 1],
                        2.0 * PI * freqs_hz[i0 + 2],
                        2.0 * PI * freqs_hz[i0 + 3],
                    ];
                    let ims = [
                        y_red_grid[i0][(a_slot, a_slot)].im,
                        y_red_grid[i0 + 1][(a_slot, a_slot)].im,
                        y_red_grid[i0 + 2][(a_slot, a_slot)].im,
                        y_red_grid[i0 + 3][(a_slot, a_slot)].im,
                    ];
                    let res = [
                        y_red_grid[i0][(a_slot, a_slot)].re,
                        y_red_grid[i0 + 1][(a_slot, a_slot)].re,
                        y_red_grid[i0 + 2][(a_slot, a_slot)].re,
                        y_red_grid[i0 + 3][(a_slot, a_slot)].re,
                    ];
                    // Guard against nearby pole distorting cubic stencil
                    let lin_slope = (y1.im - y0.im) / (w1 - w0);
                    let (_, d0) = eval_lagrange4_and_deriv(ws, ims, wm);
                    if d0 > 0.2 * lin_slope && d0 < 5.0 * lin_slope {
                        for _ in 0..8 {
                            let (fv, dv) = eval_lagrange4_and_deriv(ws, ims, wm);
                            if dv <= 0.0 {
                                break;
                            }
                            let next = wm - fv / dv;
                            if next < w0 || next > w1 {
                                break;
                            }
                            wm = next;
                        }
                        let (_, slope_m) = eval_lagrange4_and_deriv(ws, ims, wm);
                        let (re_m, _) = eval_lagrange4_and_deriv(ws, res, wm);
                        (re_m, slope_m)
                    } else {
                        let t = (wm - w0) / (w1 - w0);
                        (y0.re + t * (y1.re - y0.re), lin_slope)
                    }
                } else {
                    let t = (wm - w0) / (w1 - w0);
                    (y0.re + t * (y1.re - y0.re), (y1.im - y0.im) / (w1 - w0))
                };

                if d_im_dw > 0.0 {
                    let c_sigma_f = 0.5 * d_im_dw;
                    let f_m_hz = wm / (2.0 * PI);
                    let l_eff_h = 1.0 / (wm * wm * c_sigma_f);
                    let re_pos = re_y.max(0.0);
                    let kappa_rad_s = re_pos / c_sigma_f;
                    let kappa_hz = kappa_rad_s / (2.0 * PI);
                    let q_m = if kappa_hz > 1.0e-15 * f_m_hz { f_m_hz / kappa_hz } else { f64::INFINITY };
                    let t1_purcell_s = if kappa_rad_s > 1.0e-30 { 1.0 / kappa_rad_s } else { f64::INFINITY };
                    push_dedup(
                        &mut modes,
                        BbqMode { port: port_name.clone(), f_m_hz, c_sigma_f, l_eff_h, kappa_hz, q_m, t1_purcell_s },
                    );
                }
            }
        }

        // Fine local bracket scan around QEP eigenfrequency hints to catch narrow dispersive pole-zero pairs
        if let Some(eval_fn) = exact_rom_eval {
            let f_min_grid = freqs_hz[0] * 0.98;
            let f_max_grid = freqs_hz[n_pts - 1] * 1.02;
            for &f_hint in qep_hints_hz {
                if !(f_hint.is_finite() && f_hint >= f_min_grid && f_hint <= f_max_grid && f_hint > 0.0) {
                    continue;
                }
                if modes.iter().any(|m| m.port == *port_name && (m.f_m_hz - f_hint).abs() <= 5.0e-3 * f_hint) {
                    continue;
                }
                let n_local = 256usize;
                let w_lo = 2.0 * PI * f_hint * 0.99;
                let w_hi = 2.0 * PI * f_hint * 1.01;
                let mut prev_w = w_lo;
                let mut prev_im = eval_fn(w_lo).map(|ym| ym[(a_slot, a_slot)].im).unwrap_or(f64::NAN);
                for s in 1..=n_local {
                    let cur_w = w_lo + (w_hi - w_lo) * (s as f64) / (n_local as f64);
                    let cur_im = eval_fn(cur_w).map(|ym| ym[(a_slot, a_slot)].im).unwrap_or(f64::NAN);
                    if prev_im.is_finite() && cur_im.is_finite() && prev_im <= 0.0 && cur_im > 0.0 {
                        if let Some(m) = refine_rom_bracket(eval_fn, a_slot, port_name, prev_w, cur_w, prev_im, cur_im)
                        {
                            push_dedup(&mut modes, m);
                        }
                    }
                    prev_w = cur_w;
                    prev_im = cur_im;
                }
            }
        }
    }

    modes.sort_by(|a, b| a.f_m_hz.partial_cmp(&b.f_m_hz).unwrap());
    modes
}

// ------------------------------------------------------------------ 3. Taubin Algebraic Circle Fit + Phase Fit (Probst et al. 2015)

/// Fits a circle `(x_c, y_c, radius)` to 2D points `(x_k, y_k)` via centered `3 x 3` Taubin (1991) algebraic fit.
fn taubin_circle_fit(pts: &[[f64; 2]]) -> Result<(f64, f64, f64), Error> {
    let n = pts.len();
    if n < 3 {
        return Err(Error::invalid("Taubin S21 circle fit requires at least 3 frequency points"));
    }
    let inv_n = 1.0 / (n as f64);
    let mx: f64 = pts.iter().map(|p| p[0]).sum::<f64>() * inv_n;
    let my: f64 = pts.iter().map(|p| p[1]).sum::<f64>() * inv_n;

    let mut var_sum = 0.0;
    for p in pts {
        let dx = p[0] - mx;
        let dy = p[1] - my;
        var_sum += dx * dx + dy * dy;
    }
    let scale = (var_sum * inv_n).sqrt();
    if !(scale > 1.0e-15 && scale.is_finite()) {
        return Err(Error::invalid(
            "S21 trace has zero variance (flat off-resonant trace); cannot fit resonance circle",
        ));
    }

    // Centered and unit-RMS-scaled coordinates (u_k, v_k): mean(u) = mean(v) = 0, mean(z) = 1 where z_k = u_k^2 + v_k^2
    let (mut m_zz, mut m_zu, mut m_zv, mut m_uu, mut m_uv, mut m_vv) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut mean_z = 0.0;
    for p in pts {
        let u = (p[0] - mx) / scale;
        let v = (p[1] - my) / scale;
        let z = u * u + v * v;
        mean_z += z;
        m_zz += z * z;
        m_zu += z * u;
        m_zv += z * v;
        m_uu += u * u;
        m_uv += u * v;
        m_vv += v * v;
    }
    mean_z *= inv_n;
    m_zz = m_zz * inv_n - mean_z * mean_z;
    m_zu *= inv_n;
    m_zv *= inv_n;
    m_uu *= inv_n;
    m_uv *= inv_n;
    m_vv *= inv_n;

    // 3x3 Schur-complemented Taubin pencil: M_tilde a = eta * diag(4 * mean_z, 1, 1) a
    let m3 = [m_zz, m_zu, m_zv, m_zu, m_uu, m_uv, m_zv, m_uv, m_vv];
    let b3 = [4.0 * mean_z.max(1e-12), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let (evals, evecs) = dense_gen_eigh(&m3, &b3, 3).map_err(Error::from_kernel)?;

    // Pick eigenvector associated with smallest non-negative eigenvalue
    let k_best = (0..3).find(|&k| evals[k] >= -1.0e-12).unwrap_or(0);
    let a_vec = &evecs[k_best];
    let (a, b, c) = (a_vec[0], a_vec[1], a_vec[2]);
    if a.abs() <= 1.0e-12 * (b.hypot(c)).max(1.0) {
        return Err(Error::invalid("S21 points are nearly collinear; cannot fit resonance circle"));
    }
    let d = -a * mean_z;
    let uc = -b / (2.0 * a);
    let vc = -c / (2.0 * a);
    let r_sq = uc * uc + vc * vc - d / a;
    if !(r_sq > 0.0 && r_sq.is_finite()) {
        return Err(Error::invalid("Taubin circle fit yielded non-positive radius"));
    }
    let xc = mx + scale * uc;
    let yc = my + scale * vc;
    let radius = scale * r_sq.sqrt();
    Ok((xc, yc, radius))
}

#[inline]
fn wrap_angle(mut a: f64) -> f64 {
    while a <= -PI {
        a += 2.0 * PI;
    }
    while a > PI {
        a -= 2.0 * PI;
    }
    a
}

/// Solves a `3 x 3` real linear system `A x = b` via Gaussian elimination with partial pivoting.
fn solve_3x3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for k in 0..3 {
        let mut piv = k;
        let mut best = a[k][k].abs();
        for (i, row) in a.iter().enumerate().skip(k + 1) {
            if row[k].abs() > best {
                best = row[k].abs();
                piv = i;
            }
        }
        if best <= 1.0e-30 {
            return None;
        }
        if piv != k {
            a.swap(k, piv);
            b.swap(k, piv);
        }
        for i in (k + 1)..3 {
            let f = a[i][k] / a[k][k];
            for j in k..3 {
                a[i][j] -= f * a[k][j];
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut s = b[i];
        for j in (i + 1)..3 {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

/// Fits `S_{out,in}(f)` via Probst et al. (2015) Taubin circle fit + Gauss-Newton arctan phase fit.
fn fit_s21_circle(
    ports: [String; 2],
    branch: ResonatorBranch,
    freqs_hz: &[f64],
    trace: &[[f64; 2]],
) -> Result<S21CircleFit, Error> {
    if freqs_hz.len() != trace.len() || freqs_hz.len() < 4 {
        return Err(Error::invalid("S21 circle fit requires at least 4 frequency points"));
    }
    let n = freqs_hz.len();
    let (xc, yc, radius) = taubin_circle_fit(trace)?;

    // Unwrap phase around circle center (xc, yc)
    let mut thetas = Vec::with_capacity(n);
    let mut prev_raw = 0.0;
    let mut acc = 0.0;
    for (k, p) in trace.iter().enumerate() {
        let raw = (p[1] - yc).atan2(p[0] - xc);
        if k == 0 {
            acc = raw;
        } else {
            acc += wrap_angle(raw - prev_raw);
        }
        prev_raw = raw;
        thetas.push(acc);
    }

    // Determine winding direction s_dir in {+1, -1} of theta(f) = theta_0 + 2 * s_dir * atan(2 * Q_l * (f - f_r) / f_r)
    let total_delta = thetas[n - 1] - thetas[0];
    let s_dir = if total_delta >= 0.0 { 1.0 } else { -1.0 };

    // Initial seed from maximum phase slope |d(theta)/df|
    let mut best_idx = n / 2;
    let mut max_slope = 0.0_f64;
    for k in 0..(n - 1) {
        let df = (freqs_hz[k + 1] - freqs_hz[k]).max(1.0e-12);
        let sl = ((thetas[k + 1] - thetas[k]) / df).abs();
        if sl > max_slope {
            max_slope = sl;
            best_idx = k;
        }
    }
    let f_min = freqs_hz[0];
    let f_max = freqs_hz[n - 1];
    let span = (f_max - f_min).max(1.0);
    let mut f_r =
        if best_idx + 1 < n { 0.5 * (freqs_hz[best_idx] + freqs_hz[best_idx + 1]) } else { freqs_hz[best_idx] };
    let mut q_l = (0.25 * f_r * max_slope).clamp(1.0, 1.0e10);
    let mut theta_0 = if best_idx + 1 < n { 0.5 * (thetas[best_idx] + thetas[best_idx + 1]) } else { thetas[best_idx] };

    let cost_of = |fr: f64, ql: f64, th0: f64| -> f64 {
        let mut s = 0.0;
        for k in 0..n {
            let u = 2.0 * ql * (freqs_hz[k] - fr) / fr;
            let pred = th0 + 2.0 * s_dir * u.atan();
            let r = wrap_angle(thetas[k] - pred);
            s += r * r;
        }
        s
    };

    // Levenberg-Marquardt refinement on (f_r_scaled = f_r / f_scale, log_ql = ln(q_l), theta_0)
    let f_scale = f_r.max(1.0);
    let mut lam = 1.0e-4_f64;
    let mut cur_cost = cost_of(f_r, q_l, theta_0);

    for _iter in 0..40 {
        let mut jtj = [[0.0; 3]; 3];
        let mut jtr = [0.0; 3];
        for k in 0..n {
            let fk = freqs_hz[k];
            let u = 2.0 * q_l * (fk - f_r) / f_r;
            let pred = theta_0 + 2.0 * s_dir * u.atan();
            let res = wrap_angle(thetas[k] - pred);
            let inv_1_u2 = 1.0 / (1.0 + u * u);

            // Derivatives of pred w.r.t. (f_r / f_scale, ln(Q_l), theta_0)
            let du_dfr_scaled = -2.0 * q_l * fk / (f_r * f_r) * f_scale;
            let du_dlql = u; // q_l * du/dq_l == u
            let g = [2.0 * s_dir * inv_1_u2 * du_dfr_scaled, 2.0 * s_dir * inv_1_u2 * du_dlql, 1.0];
            for ia in 0..3 {
                jtr[ia] += g[ia] * res;
                for ib in 0..3 {
                    jtj[ia][ib] += g[ia] * g[ib];
                }
            }
        }
        let mut a_lm = jtj;
        for d in 0..3 {
            a_lm[d][d] += lam * (jtj[d][d].max(1.0e-8));
        }
        let Some(step) = solve_3x3(a_lm, jtr) else {
            lam *= 10.0;
            continue;
        };
        let cand_fr = (f_r + step[0] * f_scale).clamp(0.5 * f_min, 1.5 * f_max);
        let cand_ql = (q_l * step[1].clamp(-2.0, 2.0).exp()).clamp(0.5, 1.0e12);
        let cand_th0 = theta_0 + step[2];
        let cand_cost = cost_of(cand_fr, cand_ql, cand_th0);

        if cand_cost <= cur_cost {
            let rel_step = (step[0] * f_scale / span).abs().max(step[1].abs()).max(step[2].abs());
            f_r = cand_fr;
            q_l = cand_ql;
            theta_0 = cand_th0;
            cur_cost = cand_cost;
            lam = (lam * 0.3).max(1.0e-12);
            if rel_step < 1.0e-12 {
                break;
            }
        } else {
            lam = (lam * 5.0).min(1.0e8);
            if f_r <= 0.0 {
                break;
            }
        }
    }

    let z_c = C64::new(xc, yc);
    let e_th0 = C64::new(theta_0.cos(), theta_0.sin());
    let z_res = z_c + e_th0.scale(radius);
    let z_off = z_c - e_th0.scale(radius);

    let (q_c, q_i, phi_0_rad) = match branch {
        ResonatorBranch::Hanger => {
            let a_mag = z_off.norm().max(1.0e-15);
            let z_c_norm = z_c / z_off;
            let diff = C64::ONE - z_c_norm; // (Q_l / (2 |Q_c|)) * exp(i phi_0)
            let phi_0 = diff.im.atan2(diff.re);
            let qc_abs = (q_l * a_mag) / (2.0 * radius).max(1.0e-30);
            let cos_phi = phi_0.cos();
            let qc_real = if cos_phi.abs() > 1.0e-6 { qc_abs / cos_phi } else { qc_abs };
            let inv_qi = (1.0 / q_l) - (cos_phi / qc_abs);
            let qi = if inv_qi > 1.0e-15 / q_l { 1.0 / inv_qi } else { f64::INFINITY };
            (qc_real, qi, phi_0)
        }
        ResonatorBranch::Reflection => {
            let a_mag = z_off.norm().max(1.0e-15);
            let z_c_norm = z_c / z_off;
            let diff = C64::ONE - z_c_norm; // (Q_l / |Q_c|) * exp(i phi_0)
            let phi_0 = diff.im.atan2(diff.re);
            let qc_abs = (q_l * a_mag) / radius.max(1.0e-30);
            let cos_phi = phi_0.cos();
            let qc_real = if cos_phi.abs() > 1.0e-6 { qc_abs / cos_phi } else { qc_abs };
            let inv_qi = (1.0 / q_l) - (cos_phi / qc_abs);
            let qi = if inv_qi > 1.0e-15 / q_l { 1.0 / inv_qi } else { f64::INFINITY };
            (qc_real, qi, phi_0)
        }
        ResonatorBranch::Transmission => {
            let phi_0 = z_res.im.atan2(z_res.re);
            let s21_peak = (2.0 * radius).max(1.0e-30);
            let qc = q_l / s21_peak;
            let inv_qi = (1.0 / q_l) - (1.0 / qc);
            let qi = if inv_qi > 1.0e-15 / q_l { 1.0 / inv_qi } else { f64::INFINITY };
            (qc, qi, phi_0)
        }
    };

    let kappa_hz = f_r / q_l;
    Ok(S21CircleFit {
        ports,
        resonator_branch: branch.as_str().into(),
        f_r_hz: f_r,
        q_l,
        q_c,
        q_i,
        kappa_hz,
        phi_0_rad,
        circle_center: [xc, yc],
        circle_radius: radius,
    })
}

/// Computes the 2-port (or 1-port) scattering matrix across `sel_ports` when all remaining ports `q` are
/// either left open-circuited (`lumped_jj` with `z0_ohm == None` and not in `terminations_ohm`) or terminated
/// in `terminations_ohm[q]` / `z0_ohm[q]`.
fn extract_terminated_sub_s(
    ports: &[String],
    rom_ports: Option<&[RomPort]>,
    z_full: &Dense<C64>,
    z0_row: &[f64],
    terminations_ohm: &BTreeMap<String, f64>,
    sel: [usize; 2],
) -> C64 {
    let p_tot = ports.len();
    let mut sub_ports = vec![sel[0]];
    if sel[1] != sel[0] {
        sub_ports.push(sel[1]);
    }

    // Collect external non-selected ports that have a finite resistive termination
    let mut load_idx = Vec::new();
    let mut load_z = Vec::new();
    for q in 0..p_tot {
        if sub_ports.contains(&q) {
            continue;
        }
        if let Some(&r_override) = terminations_ohm.get(&ports[q]) {
            load_idx.push(q);
            load_z.push(r_override);
        } else if let Some(rp_list) = rom_ports {
            // In CircuitRom, nonlinear lumped_jj ports without z0_ohm are internal (open-circuited externally)
            if let Some(z0) = rp_list[q].z0_ohm.filter(|&z| z.is_finite() && z > 0.0) {
                load_idx.push(q);
                load_z.push(z0);
            }
        } else if let Some(&z0) = z0_row.get(q).filter(|&&z| z.is_finite() && z > 0.0) {
            load_idx.push(q);
            load_z.push(z0);
        }
    }

    let ns = sub_ports.len();
    let nl = load_idx.len();
    let z_ss = if nl == 0 {
        Dense::from_fn(ns, ns, |i, j| z_full[(sub_ports[i], sub_ports[j])])
    } else {
        let z_sa = Dense::from_fn(ns, ns, |i, j| z_full[(sub_ports[i], sub_ports[j])]);
        let z_sl = Dense::from_fn(ns, nl, |i, j| z_full[(sub_ports[i], load_idx[j])]);
        let z_ls = Dense::from_fn(nl, ns, |i, j| z_full[(load_idx[i], sub_ports[j])]);
        let mut z_ll = Dense::from_fn(nl, nl, |i, j| z_full[(load_idx[i], load_idx[j])]);
        for (t, &zl) in load_z.iter().enumerate() {
            z_ll[(t, t)] += C64::new(zl, 0.0);
        }
        let z_ll_inv = invert_complex_matrix(&z_ll);
        let corr = z_sl.matmul(&z_ll_inv.matmul(&z_ls));
        Dense::from_fn(ns, ns, |i, j| z_sa[(i, j)] - corr[(i, j)])
    };

    let z0_sub: Vec<f64> = sub_ports
        .iter()
        .map(|&pi| terminations_ohm.get(&ports[pi]).copied().or_else(|| z0_row.get(pi).copied()).unwrap_or(50.0))
        .collect();
    let s_sub = z_to_s_matrix(&z_ss, &z0_sub);
    if sel[0] == sel[1] {
        s_sub[(0, 0)]
    } else {
        // S_{out, in}: sel[0] is input port (col 0), sel[1] is output port (row 1)
        s_sub[(1, 0)]
    }
}

// ------------------------------------------------------------------ 4. Touchstone v1.1 (.sNp) Formatter

fn format_touchstone_pair(z: C64, fmt: TouchstoneFormat) -> (f64, f64) {
    match fmt {
        TouchstoneFormat::Ri => (z.re, z.im),
        TouchstoneFormat::Ma => (z.norm(), z.im.atan2(z.re).to_degrees()),
        TouchstoneFormat::Db => (20.0 * z.norm().max(1.0e-300).log10(), z.im.atan2(z.re).to_degrees()),
    }
}

fn format_touchstone_snp(
    ports: &[String],
    freqs_hz: &[f64],
    s_mats: &[Dense<C64>],
    r_ref_ohm: f64,
    fmt: TouchstoneFormat,
) -> String {
    let n_ports = ports.len();
    let mut out = String::new();
    let _ = writeln!(out, "! Touchstone v1.1 ({}-port) generated by nyx-analysis", n_ports);
    let port_labels: Vec<String> = ports.iter().enumerate().map(|(i, p)| format!("{}={p}", i + 1)).collect();
    let _ = writeln!(out, "! Ports: {}", port_labels.join(", "));
    let _ = writeln!(out, "# Hz S {} R {:.6}", fmt.as_str(), r_ref_ohm);

    for (k, &f_hz) in freqs_hz.iter().enumerate() {
        let s = &s_mats[k];
        match n_ports {
            0 => {}
            1 => {
                let (a, b) = format_touchstone_pair(s[(0, 0)], fmt);
                let _ = writeln!(out, "{f_hz:.9e} {a:.9e} {b:.9e}");
            }
            2 => {
                // Touchstone v1.1 2-port order: S11 S21 S12 S22
                let (a11, b11) = format_touchstone_pair(s[(0, 0)], fmt);
                let (a21, b21) = format_touchstone_pair(s[(1, 0)], fmt);
                let (a12, b12) = format_touchstone_pair(s[(0, 1)], fmt);
                let (a22, b22) = format_touchstone_pair(s[(1, 1)], fmt);
                let _ = writeln!(
                    out,
                    "{f_hz:.9e} {a11:.9e} {b11:.9e} {a21:.9e} {b21:.9e} {a12:.9e} {b12:.9e} {a22:.9e} {b22:.9e}"
                );
            }
            _ => {
                // Touchstone v1.1 N >= 3 order: row-major (row i = 0..N, col j = 0..N), max 4 pairs per line
                for i in 0..n_ports {
                    for chunk_start in (0..n_ports).step_by(4) {
                        if i == 0 && chunk_start == 0 {
                            let _ = write!(out, "{f_hz:.9e}");
                        } else {
                            let _ = write!(out, "               ");
                        }
                        let chunk_end = (chunk_start + 4).min(n_ports);
                        for j in chunk_start..chunk_end {
                            let (a, b) = format_touchstone_pair(s[(i, j)], fmt);
                            let _ = write!(out, " {a:.9e} {b:.9e}");
                        }
                        let _ = writeln!(out);
                    }
                }
            }
        }
    }
    out
}

// ------------------------------------------------------------------ Public Entry Points & Analysis Trait

/// Runs [`NetworkAnalysis`] with default parameters.
pub fn network(result: &SolveResult) -> Result<NetworkReport, Error> {
    network_with_params(result, &Value::Null)
}

/// Runs [`NetworkAnalysis`] with optional JSON parameters.
pub fn network_with_params(result: &SolveResult, params: &Value) -> Result<NetworkReport, Error> {
    let t0 = Instant::now();
    let rom_opt = CircuitRom::of(result);
    let driven_opt = DrivenPortData::of(result);
    if rom_opt.is_none() && driven_opt.is_none() {
        return Err(Error::invalid(
            "the result carries neither a CircuitRom nor DrivenPortData (run an ems cavity or driven study)",
        ));
    }

    let net_params = NetworkParams::from_json(params)?;

    // Resolve CircuitRom ports + QEP if present
    let mut resolved_rom_ports: Option<Vec<RomPort>> = None;
    let mut c_internal_rom: Vec<Vec<f64>> = Vec::new();
    let mut k_tot_rom: Vec<Vec<f64>> = Vec::new();
    let mut m_tot_rom: Vec<Vec<f64>> = Vec::new();
    let mut qep_modes: Vec<QepMode> = Vec::new();

    if let Some(ref rom) = rom_opt {
        let (ports_ov, _) = resolve_ports_with_overrides(&rom.ports, params)?;
        for k in net_params.terminations_ohm.keys() {
            if driven_opt.is_none() && !ports_ov.iter().any(|p| &p.name == k) {
                return Err(Error::invalid(format!(
                    "terminations_ohm port '{k}' does not match any port in {:?}",
                    ports_ov.iter().map(|p| p.name.as_str()).collect::<Vec<_>>()
                )));
            }
        }
        let (c_int, c_tot, _) = assemble_rom_damping(rom, &ports_ov, Some(&net_params.terminations_ohm));
        qep_modes = solve_rom_qep(rom, &ports_ov, &c_tot)?;
        let (kt, mt) = assemble_rom_km(rom, &ports_ov);
        k_tot_rom = kt;
        m_tot_rom = mt;
        c_internal_rom = c_int;
        resolved_rom_ports = Some(ports_ov);
    }

    // Choose multi-port frequency-domain source: DrivenPortData (when present & no custom freqs_hz override on ROM) or CircuitRom synthesis
    let use_rom_synthesis = rom_opt.is_some() && (driven_opt.is_none() || net_params.freqs_hz.is_some());

    let (ports, freqs_hz, z0_ohm, z_mats, y_mats, s_mats) = if !use_rom_synthesis {
        let drv = driven_opt.as_ref().unwrap();
        for k in net_params.terminations_ohm.keys() {
            if !drv.ports.iter().any(|p| p == k) {
                return Err(Error::invalid(format!(
                    "terminations_ohm port '{k}' does not match any port in {:?}",
                    drv.ports
                )));
            }
        }
        let n_f = drv.freqs_hz.len();
        let p = drv.ports.len();
        let mut z_list = Vec::with_capacity(n_f);
        let mut y_list = Vec::with_capacity(n_f);
        let mut s_list = Vec::with_capacity(n_f);
        let mut z0_out = Vec::with_capacity(n_f);

        for k in 0..n_f {
            let z_dense = pairs_to_dense(&drv.z_matrix[k]);
            let y_dense = pairs_to_dense(&drv.y_matrix[k]);
            let z0_src = drv.z0_ohm.get(k).or_else(|| drv.z0_ohm.first());
            let mut z0_row_pairs = Vec::with_capacity(p);
            let mut z0_re = Vec::with_capacity(p);
            let mut has_term_override = false;
            for (pi, p_name) in drv.ports.iter().enumerate() {
                if let Some(&r_ov) = net_params.terminations_ohm.get(p_name) {
                    z0_row_pairs.push([r_ov, 0.0]);
                    z0_re.push(r_ov);
                    has_term_override = true;
                } else {
                    let pair = z0_src.and_then(|r| r.get(pi)).copied().unwrap_or([50.0, 0.0]);
                    z0_row_pairs.push(pair);
                    z0_re.push(pair[0].max(1.0e-12));
                }
            }
            let s_dense =
                if has_term_override { z_to_s_matrix(&z_dense, &z0_re) } else { pairs_to_dense(&drv.s_matrix[k]) };
            z0_out.push(z0_row_pairs);
            z_list.push(z_dense);
            y_list.push(y_dense);
            s_list.push(s_dense);
        }
        (drv.ports.clone(), drv.freqs_hz.clone(), z0_out, z_list, y_list, s_list)
    } else {
        let rom = rom_opt.as_ref().unwrap();
        let rports = resolved_rom_ports.as_ref().unwrap();
        let port_names: Vec<String> = rports.iter().map(|p| p.name.clone()).collect();
        let p = port_names.len();

        // Synthesize on explicit `freqs_hz`, or auto-grid around QEP modes if BBQ / S21 / Touchstone was requested
        let grid: Vec<f64> = if let Some(ref f_user) = net_params.freqs_hz {
            f_user.clone()
        } else if net_params.active_ports.is_some()
            || net_params.s21_ports.is_some()
            || net_params.touchstone_format.is_some()
        {
            if !qep_modes.is_empty() {
                let f_lo = 0.6 * qep_modes.first().unwrap().f_re_hz;
                let f_hi = 1.3 * qep_modes.last().unwrap().f_re_hz;
                let n_auto = 251usize;
                (0..n_auto).map(|k| f_lo + (f_hi - f_lo) * (k as f64) / ((n_auto - 1) as f64)).collect()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };

        let mut z0_re = Vec::with_capacity(p);
        let mut z0_pairs_single = Vec::with_capacity(p);
        for rp in rports {
            let z0 = net_params.terminations_ohm.get(&rp.name).copied().or(rp.z0_ohm).unwrap_or(50.0);
            z0_re.push(z0);
            z0_pairs_single.push([z0, 0.0]);
        }

        let mut z_list = Vec::with_capacity(grid.len());
        let mut y_list = Vec::with_capacity(grid.len());
        let mut s_list = Vec::with_capacity(grid.len());
        let mut z0_out = Vec::with_capacity(grid.len());

        for &f_hz in &grid {
            let omega = 2.0 * PI * f_hz;
            let (z_m, y_m) = eval_rom_zy_at_omega(rom, &k_tot_rom, &m_tot_rom, &c_internal_rom, omega)
                .ok_or_else(|| Error::invalid(format!("CircuitRom Z/Y synthesis is singular at {f_hz} Hz")))?;
            let s_m = z_to_s_matrix(&z_m, &z0_re);
            z0_out.push(z0_pairs_single.clone());
            z_list.push(z_m);
            y_list.push(y_m);
            s_list.push(s_m);
        }
        (port_names, grid, z0_out, z_list, y_list, s_list)
    };

    // Validate active_ports and s21_ports against `ports`
    let port_index = |name: &str, ctx: &str| -> Result<usize, Error> {
        ports
            .iter()
            .position(|p| p == name)
            .ok_or_else(|| Error::invalid(format!("{ctx} port '{name}' does not match any port in {ports:?}")))
    };

    let active_idx: Vec<usize> = if let Some(ref ap) = net_params.active_ports {
        let mut seen = BTreeSet::new();
        let mut out = Vec::with_capacity(ap.len());
        for name in ap {
            let idx = port_index(name, "active_ports")?;
            if seen.insert(idx) {
                out.push(idx);
            }
        }
        out
    } else if use_rom_synthesis && !freqs_hz.is_empty() {
        resolved_rom_ports
            .as_ref()
            .unwrap()
            .iter()
            .enumerate()
            .filter(|(_, p)| is_nonlinear_port(p))
            .map(|(i, _)| i)
            .collect()
    } else if !use_rom_synthesis && !net_params.terminations_ohm.is_empty() {
        (0..ports.len()).filter(|&i| !net_params.terminations_ohm.contains_key(&ports[i])).collect()
    } else if !use_rom_synthesis && ports.len() == 1 && net_params.s21_ports.is_none() {
        vec![0]
    } else {
        Vec::new()
    };

    let active_set: BTreeSet<usize> = active_idx.iter().copied().collect();
    // Terminated ports T for Schur BBQ: exclude internal open-circuited `lumped_jj` ports without z0_ohm/termination
    let term_idx: Vec<usize> = (0..ports.len())
        .filter(|&i| {
            if active_set.contains(&i) {
                return false;
            }
            if net_params.terminations_ohm.contains_key(&ports[i]) {
                return true;
            }
            if let Some(ref rports) = resolved_rom_ports {
                if use_rom_synthesis {
                    return rports[i].z0_ohm.is_some_and(|z| z.is_finite() && z > 0.0);
                }
            }
            true
        })
        .collect();

    let z0_grid_re: Vec<Vec<f64>> =
        z0_ohm.iter().map(|row| row.iter().map(|pair| pair[0].max(1.0e-12)).collect()).collect();

    let bbq_modes = if !active_idx.is_empty() && freqs_hz.len() >= 2 {
        if use_rom_synthesis {
            let rom_ref = rom_opt.as_ref().unwrap();
            let rports_ref = resolved_rom_ports.as_ref().unwrap();
            let qep_hints: Vec<f64> = qep_modes.iter().map(|m| m.f_re_hz).collect();
            let z_loads_rom: Vec<f64> = term_idx
                .iter()
                .map(|&ti| {
                    net_params.terminations_ohm.get(&ports[ti]).copied().or(rports_ref[ti].z0_ohm).unwrap_or(50.0)
                })
                .collect();
            let eval_exact = |omega: f64| -> Option<Dense<C64>> {
                let (zm, ym) = eval_rom_zy_at_omega(rom_ref, &k_tot_rom, &m_tot_rom, &c_internal_rom, omega)?;
                Some(schur_reduced_admittance(&ym, &zm, false, &active_idx, &term_idx, &z_loads_rom))
            };
            extract_bbq_modes(
                &ports,
                &freqs_hz,
                &y_mats,
                &z_mats,
                &z0_grid_re,
                &active_idx,
                &term_idx,
                &net_params.terminations_ohm,
                &qep_hints,
                Some(&eval_exact),
            )
        } else {
            extract_bbq_modes(
                &ports,
                &freqs_hz,
                &y_mats,
                &z_mats,
                &z0_grid_re,
                &active_idx,
                &term_idx,
                &net_params.terminations_ohm,
                &[],
                None,
            )
        }
    } else {
        Vec::new()
    };

    // Taubin S21 circle fit
    let resolved_s21_ports: Option<[usize; 2]> = if let Some([ref p0, ref p1]) = net_params.s21_ports {
        Some([port_index(p0, "s21_ports[0]")?, port_index(p1, "s21_ports[1]")?])
    } else if net_params.has_explicit_branch && !freqs_hz.is_empty() {
        match (net_params.resonator_branch, ports.len()) {
            (ResonatorBranch::Reflection, p) if p >= 1 => Some([0, 0]),
            (_, p) if p >= 2 => Some([0, 1]),
            _ => None,
        }
    } else {
        None
    };

    let s21_fit = if let Some([idx_in, idx_out]) = resolved_s21_ports {
        if freqs_hz.len() < 4 {
            return Err(Error::invalid("s21_ports circle fit requires at least 4 frequency points"));
        }
        let rom_slice = if use_rom_synthesis { resolved_rom_ports.as_deref() } else { None };
        let trace: Vec<[f64; 2]> = (0..freqs_hz.len())
            .map(|k| {
                let sk = extract_terminated_sub_s(
                    &ports,
                    rom_slice,
                    &z_mats[k],
                    &z0_grid_re[k],
                    &net_params.terminations_ohm,
                    [idx_in, idx_out],
                );
                [sk.re, sk.im]
            })
            .collect();
        Some(fit_s21_circle(
            [ports[idx_in].clone(), ports[idx_out].clone()],
            net_params.resonator_branch,
            &freqs_hz,
            &trace,
        )?)
    } else {
        None
    };

    // Touchstone v1.1 (.sNp) formatting
    let touchstone_snp = if !freqs_hz.is_empty() && !ports.is_empty() {
        let fmt = net_params.touchstone_format.unwrap_or(TouchstoneFormat::Ri);
        let first_z0 = z0_grid_re.first().and_then(|r| r.first()).copied().unwrap_or(50.0);
        let uniform_z0 =
            z0_grid_re.iter().all(|row| row.iter().all(|&z| (z - first_z0).abs() <= 1.0e-9 * first_z0.max(1.0)));
        let r_ref = if uniform_z0 { first_z0 } else { 50.0 };
        let s_for_snp: Vec<Dense<C64>> = if uniform_z0 {
            s_mats.clone()
        } else {
            let ref_vec = vec![r_ref; ports.len()];
            z_mats.iter().map(|zm| z_to_s_matrix(zm, &ref_vec)).collect()
        };
        Some(format_touchstone_snp(&ports, &freqs_hz, &s_for_snp, r_ref, fmt))
    } else {
        None
    };

    let s_matrix = s_mats.iter().map(dense_to_pairs).collect();
    let z_matrix = z_mats.iter().map(dense_to_pairs).collect();
    let y_matrix = y_mats.iter().map(dense_to_pairs).collect();

    Ok(NetworkReport {
        ports,
        freqs_hz,
        z0_ohm,
        s_matrix,
        z_matrix,
        y_matrix,
        qep_modes,
        bbq_modes,
        s21_fit,
        touchstone_snp,
        elapsed_s: t0.elapsed().as_secs_f64(),
    })
}

/// Unified multi-port microwave network, QEP, BBQ, `S_21` fit, and Touchstone analysis (`id: "network"`).
pub struct NetworkAnalysis;

impl Analysis for NetworkAnalysis {
    fn id(&self) -> &'static str {
        "network"
    }

    fn name(&self) -> &'static str {
        "Multi-port microwave network, QEP, Black-Box Quantization (BBQ), S21 fit, and Touchstone"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        network_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        Ok(network(result)?.to_json())
    }

    fn run_with_params(&self, _mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        Ok(network_with_params(result, params)?.to_json())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::epr::jj_flux_scale;
    use crate::hamiltonian::hamiltonian_with_params;

    #[test]
    fn qep_schur_bbq_and_hamiltonian_purcell_t1_agree_on_transmon_resonator_network() {
        // Transmon (f_q ~ 4.8 GHz, C_q = 75 fF) + readout resonator (f_r ~ 6.8 GHz, C_r = 50 fF) + 2 feedline ports (50 Ohm each)
        let ej_hz = 15.0e9;
        let lj_h = jj_flux_scale() / ej_hz;
        let c_q = 75.0e-15;
        let c_r = 50.0e-15;
        let c_g = 3.0e-15;
        let f_r_bare = 6.8e9;
        let omega_r = 2.0 * PI * f_r_bare;
        let l_r_inv = omega_r * omega_r * (c_r + c_g);

        // Coupled to 2 feedline ports (inline) with Z0 = 50 Ohm so loaded Q_l ~ 10_000 (kappa/2pi ~ 0.68 MHz)
        let b_feed = ((omega_r * (c_r + c_g) * 50.0) / 20000.0).sqrt();

        let rom = CircuitRom {
            ports: vec![
                RomPort {
                    name: "Q1_JJ".into(),
                    kind: "lumped_jj".into(),
                    l_j_h: Some(lj_h),
                    e_j_hz: Some(ej_hz),
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: None,
                    z0_ohm: None,
                    direction: [0.0, 1.0, 0.0],
                    formulation: "rank1".into(),
                },
                RomPort {
                    name: "feed_1".into(),
                    kind: "lumped_port".into(),
                    l_j_h: None,
                    e_j_hz: None,
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: None,
                    z0_ohm: Some(50.0),
                    direction: [1.0, 0.0, 0.0],
                    formulation: "rank1".into(),
                },
                RomPort {
                    name: "feed_2".into(),
                    kind: "lumped_port".into(),
                    l_j_h: None,
                    e_j_hz: None,
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: None,
                    z0_ohm: Some(50.0),
                    direction: [1.0, 0.0, 0.0],
                    formulation: "rank1".into(),
                },
            ],
            r_static: 0,
            r_dyn: 2,
            k_r: vec![vec![0.0, 0.0], vec![0.0, MU0 * l_r_inv]],
            c_r: vec![],
            m_r: vec![vec![(c_q + c_g) / EPS0, -c_g / EPS0], vec![-c_g / EPS0, (c_r + c_g) / EPS0]],
            b_r: vec![vec![1.0, 0.0, 0.0], vec![0.0, b_feed, b_feed]],
            c_port_f: vec![],
            modes: vec![],
        };

        let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };
        assert!(NetworkAnalysis.applies_to(&res));

        // Non-uniform ascending sweep covering [4.0, 7.5] GHz with fine VNA resolution across the ~6.81 GHz resonator
        let mut freqs: Vec<f64> = (0..151).map(|k| 4.0e9 + (7.5e9 - 4.0e9) * (k as f64) / 150.0).collect();
        for k in 0..150 {
            freqs.push(6.805e9 + (6.825e9 - 6.805e9) * (k as f64) / 149.0);
        }
        freqs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        freqs.dedup_by(|a, b| (*a - *b).abs() < 1.0);
        let n_freqs = freqs.len();

        let net = network_with_params(
            &res,
            &json!({
                "active_ports": ["Q1_JJ"],
                "s21_ports": ["feed_1", "feed_2"],
                "resonator_branch": "transmission",
                "freqs_hz": freqs,
                "touchstone_format": "DB"
            }),
        )
        .unwrap();

        assert_eq!(net.qep_modes.len(), 2);
        assert_eq!(net.bbq_modes.len(), 2);

        // 1. QEP vs Schur BBQ mode frequencies & linewidths match to < 0.2%
        for m in 0..2 {
            let rel_f = (net.qep_modes[m].f_re_hz - net.bbq_modes[m].f_m_hz).abs() / net.qep_modes[m].f_re_hz;
            let rel_k = (net.qep_modes[m].kappa_hz - net.bbq_modes[m].kappa_hz).abs() / net.qep_modes[m].kappa_hz;
            assert!(rel_f < 5.0e-4, "mode {m} QEP vs BBQ freq mismatch: {rel_f:e}");
            assert!(rel_k < 2.0e-3, "mode {m} QEP vs BBQ kappa mismatch: {rel_k:e}");
        }

        // 2. S21 inline transmission circle fit recovers dressed cavity frequency & loaded Q_l within < 0.2%
        let fit = net.s21_fit.as_ref().unwrap();
        let rel_s21_f = (fit.f_r_hz - net.qep_modes[1].f_re_hz).abs() / net.qep_modes[1].f_re_hz;
        let rel_s21_q = (fit.q_l - net.qep_modes[1].q_c).abs() / net.qep_modes[1].q_c;
        assert!(rel_s21_f < 1.0e-4, "S21 fit f_r err = {rel_s21_f:e}");
        assert!(
            rel_s21_q < 2.0e-3,
            "S21 fit Q_l ({}) vs QEP Q ({}): err = {rel_s21_q:e}",
            fit.q_l,
            net.qep_modes[1].q_c
        );

        // 3. Touchstone v1.1 .s3p formatting check
        let snp = net.touchstone_snp.as_ref().unwrap();
        assert!(snp.contains("3-port"));
        assert!(snp.contains("# Hz S DB R 50.000000"));
        let data_lines: Vec<&str> = snp.lines().filter(|l| !l.starts_with('!') && !l.starts_with('#')).collect();
        // 3-port .s3p emits 3 lines per frequency point (one row per port)
        assert_eq!(data_lines.len(), n_freqs * 3);

        // 4. FR-22 cross-check: bare `kappa_cav_Hz` equals exact sum of dressed QEP linewidths (trace invariance)
        //    and agrees with dressed cavity linewidth within dispersive hybridization `(g/Delta)^2 ~ 0.7%`!
        let ham_cb =
            hamiltonian_with_params(&res, &json!({ "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [5] }))
                .unwrap();
        let sum_qep_kappa = net.qep_modes[0].kappa_hz + net.qep_modes[1].kappa_hz;
        let rel_trace = (ham_cb.kappa_cav_hz[0] - sum_qep_kappa).abs() / ham_cb.kappa_cav_hz[0];
        assert!(
            rel_trace < 1.0e-5,
            "trace invariance mismatch: bare={} vs sum_qep={}",
            ham_cb.kappa_cav_hz[0],
            sum_qep_kappa
        );
        let rel_kcav = (ham_cb.kappa_cav_hz[0] - net.qep_modes[1].kappa_hz).abs() / net.qep_modes[1].kappa_hz;
        assert!(
            rel_kcav < 1.0e-2,
            "kappa_cav_Hz ({}) vs QEP ({}) rel = {rel_kcav:e}",
            ham_cb.kappa_cav_hz[0],
            net.qep_modes[1].kappa_hz
        );

        // Note: ham_cb evaluates Purcell T1 at the anharmonic f01 (~4.79 GHz vs harmonic 5.06 GHz); check harmonic limit E_J -> large or compare with exact Schur at f01
        let ej_harm = 120.0e9;
        let ham_harm = hamiltonian_with_params(
            &res,
            &json!({ "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [5], "ej_hz_override": { "Q1_JJ": ej_harm } }),
        )
        .unwrap();
        let freqs_hi: Vec<f64> = (0..301).map(|k| 6.0e9 + (16.0e9 - 6.0e9) * (k as f64) / 300.0).collect();
        let net_harm = network_with_params(
            &res,
            &json!({
                "active_ports": ["Q1_JJ"],
                "freqs_hz": freqs_hi,
                "ej_hz_override": { "Q1_JJ": ej_harm }
            }),
        )
        .unwrap();
        let bbq_qubit = net_harm.bbq_modes.iter().max_by(|a, b| a.f_m_hz.partial_cmp(&b.f_m_hz).unwrap()).unwrap();
        let rel_t1 = (ham_harm.t1_purcell_s[0] - bbq_qubit.t1_purcell_s).abs() / bbq_qubit.t1_purcell_s;
        assert!(
            rel_t1 < 0.03,
            "Hamiltonian t1_purcell_s ({:e}) vs Schur BBQ ({:e}): rel = {rel_t1:.4}",
            ham_harm.t1_purcell_s[0],
            bbq_qubit.t1_purcell_s
        );
    }

    #[test]
    fn taubin_hanger_and_reflection_circle_fits_recover_exact_tilted_resonator_parameters() {
        // Synthetic tilted hanger notch: f_r = 6.5 GHz, Q_i = 200_000, |Q_c| = 25_000, phi_0 = 0.22 rad
        let f_r = 6.5e9_f64;
        let q_i = 200_000.0_f64;
        let qc_abs = 25_000.0_f64;
        let phi_0 = 0.22_f64;
        let inv_ql = (1.0 / q_i) + (phi_0.cos() / qc_abs);
        let q_l = 1.0 / inv_ql;

        let n_pts = 81usize;
        let span = 6.0 * (f_r / q_l);
        let freqs_hz: Vec<f64> =
            (0..n_pts).map(|k| (f_r - 0.5 * span) + span * (k as f64) / ((n_pts - 1) as f64)).collect();

        let mut s_matrix = Vec::with_capacity(n_pts);
        let mut z_matrix = Vec::with_capacity(n_pts);
        let mut y_matrix = Vec::with_capacity(n_pts);
        let mut z0_ohm = Vec::with_capacity(n_pts);

        for &f in &freqs_hz {
            let denom = C64::new(1.0, 2.0 * q_l * (f - f_r) / f_r);
            let num = C64::new(phi_0.cos(), phi_0.sin()).scale(q_l / qc_abs);
            let s21 = C64::ONE - num / denom;
            let s11 = -(num / denom);
            let mut s_mat = Dense::<C64>::zeros(2, 2);
            s_mat[(0, 0)] = s11;
            s_mat[(1, 1)] = s11;
            s_mat[(1, 0)] = s21;
            s_mat[(0, 1)] = s21;
            // Z = 50 (I + S)(I - S)^-1
            let i_minus_s = Dense::from_fn(2, 2, |i, j| (if i == j { C64::ONE } else { C64::ZERO }) - s_mat[(i, j)]);
            let i_plus_s = Dense::from_fn(2, 2, |i, j| (if i == j { C64::ONE } else { C64::ZERO }) + s_mat[(i, j)]);
            let z_mat = i_plus_s.matmul(&invert_complex_matrix(&i_minus_s));
            let z_mat = Dense::from_fn(2, 2, |i, j| z_mat[(i, j)].scale(50.0));
            let y_mat = invert_complex_matrix(&z_mat);

            s_matrix.push(dense_to_pairs(&s_mat));
            z_matrix.push(dense_to_pairs(&z_mat));
            y_matrix.push(dense_to_pairs(&y_mat));
            z0_ohm.push(vec![[50.0, 0.0], [50.0, 0.0]]);
        }

        let drv = DrivenPortData {
            freqs_hz,
            ports: vec!["in".into(), "out".into()],
            z0_ohm,
            s_matrix,
            z_matrix,
            y_matrix,
            hdm_solves: 3,
            rom_dim: 6,
            anchor_freqs_hz: vec![f_r],
            max_residual: 1e-9,
        };
        let res = SolveResult { scalars: vec![drv.scalar()], ..Default::default() };

        let rep = network_with_params(
            &res,
            &json!({
                "s21_ports": ["in", "out"],
                "resonator_branch": "hanger",
                "touchstone_format": "RI"
            }),
        )
        .unwrap();

        let fit = rep.s21_fit.unwrap();
        assert!((fit.f_r_hz - f_r).abs() / f_r < 1.0e-10, "f_r err: {}", fit.f_r_hz - f_r);
        assert!((fit.q_l - q_l).abs() / q_l < 1.0e-8, "q_l got {} want {}", fit.q_l, q_l);
        assert!((fit.q_i - q_i).abs() / q_i < 1.0e-6, "q_i got {} want {}", fit.q_i, q_i);
        assert!((fit.phi_0_rad - phi_0).abs() < 1.0e-8, "phi_0 got {} want {}", fit.phi_0_rad, phi_0);
        assert!(rep.touchstone_snp.unwrap().contains("# Hz S RI R 50.000000"));
    }
}
