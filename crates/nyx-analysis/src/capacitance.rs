// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Multi-terminal capacitance (`capacitance_matrix`) and inductance (`inductance_matrix`) analyses of `dfs` results.
//!
//! `CapacitanceMatrixAnalysis` reads `CapacitanceMatrixData` from `SolveResult::scalars` (produced by an electrostatic
//! study of `dfs` with one or more `terminal` boundary conditions) and returns the Maxwell capacitance matrix (F), the
//! mutual/SPICE capacitance matrix (F), and the charging energy matrix `E_C / h = (e^2 / 2h) C_M^(-1)` (Hz).
//!
//! When optional `params` are supplied to [`capacitance_matrix_with_mesh_and_params`] (or `Analysis::run_with_params`),
//! it also performs:
//! - **3D Matched-Asymptotic Finite-Thickness Correction (`metal_thickness_m`, FR-11 Pt 1)**: Integrates the exact
//!   multi-terminal singular edge-tube Gram matrix `\tilde{I}_{edge,ij} = (1 / R_0) \int_{r <= R_0} \eps E^(i) . E^(j) dV`
//!   around the 1D perimeter `\partial\Gamma_{SC}` of 3D zero-thickness `Op::Sheet` conductors (indexed via `SegmentBvh3`),
//!   evaluates the intrinsic gap scale `s_eff = tr(C^M) / tr(\tilde{I}_{edge})`, and applies the Schwarz–Christoffel
//!   matched-asymptotic capacitance shift `\Delta C^M(t) = t [c_0 + c_1 ln(s_eff / t)] \tilde{I}_{edge}` (asserting SPD)
//!   prior to Schur/qubit reduction;
//! - Ground-terminal projection (`"eliminate": "drop"`) or floating-reference Schur elimination
//!   (`"eliminate": "floating_schur"`, `C_ss - C_sg C_gg^-1 C_gs`), reporting `reduced_terminals`,
//!   `maxwell_reduced_F`, `mutual_reduced_F`, `ec_reduced_Hz`, `condition_number`, `null_modes`, and `has_null_mode`;
//! - Single-pad grounded (`Xmon`) and differential two-pad floating (`TransmonPocket`, Eqs. 17–18 of
//!   `arXiv:2611.22057v1`) lumped qubit extraction (`C_Sigma_F`, `E_C_Hz`, `C_eff_F`, `beta`, and transmon `ej_over_ec`,
//!   `f01_Hz`, `anharmonicity_Hz`, `n_zpf` when `E_J_Hz` is given);
//! - Quarter-wave, half-wave, and lumped readout/bus resonator ZPF voltage (`C_r_F`, `V_rms_V`) and pairwise
//!   capacitive coupling `g_Hz` and dispersive shift `chi_Hz`.
//!
//! `InductanceMatrixAnalysis` reads `InductanceMatrixData` from `SolveResult::scalars` (produced by a magnetostatic
//! study of `dfs`) and returns the multi-terminal inductance matrix `L` (H, or H/m for 2D `cross_section` models) and
//! the active magnetostatic formulation (`"cross_section"`, `"planar_stream"`, or `"mixed_3d"`).
#![allow(clippy::needless_range_loop)]

use crate::epr::{EPS0, E_CHARGE, PLANCK_H};
use crate::tls::{detect_sheet_perimeter_3d, pt3, subdivided_tet_quad, EFieldAccessor};
use crate::Analysis;
use femcore::fem::{geometry, quad_rule, shape};
use femcore::linalg::driver::dense_eigh;
use nyx_physics::{discretise, CapacitanceMatrixData, Error, InductanceMatrixData, Mesh, SolveResult, TlsMeta};
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

/// Relative eigenvalue floor below which a symmetric capacitance mode is classified as a gauge/floating null mode.
const NULL_MODE_REL_TOL: f64 = 1.0e-10;
/// Finite ceiling used for `condition_number` when an exact zero/floating null mode is present (keeping JSON finite).
const COND_CAP: f64 = 1.0e16;

/// Extracts the typed [`CapacitanceMatrixData`] payload from `result`.
pub fn capacitance_matrix(result: &SolveResult) -> Result<CapacitanceMatrixData, Error> {
    CapacitanceMatrixData::of(result).ok_or_else(|| {
        Error::invalid(
            "the result carries no capacitance matrix (run an electrostatic study of dfs with terminal conditions)",
        )
    })
}

/// True when `result` carries a multi-terminal capacitance matrix payload.
pub fn capacitance_matrix_applies(result: &SolveResult) -> bool {
    CapacitanceMatrixData::of(result).is_some()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EliminateMode {
    Drop,
    FloatingSchur,
}

impl EliminateMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::Drop => "drop",
            Self::FloatingSchur => "floating_schur",
        }
    }
}

struct ComputedQubit {
    name: String,
    pads: Vec<String>,
    c_sigma_f: f64,
    e_c_hz: f64,
    c_eff_f: BTreeMap<String, f64>,
    beta: BTreeMap<String, f64>,
    e_j_hz: Option<f64>,
    f01_hz: Option<f64>,
    n_zpf: Option<f64>,
    json: Value,
}

struct ComputedResonator {
    name: String,
    terminal: String,
    f_r_hz: f64,
    v_rms_v: f64,
    json: Value,
}

fn mutual_from_maxwell(maxwell: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = maxwell.len();
    let mut mutual = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut row_sum = 0.0;
        for (k, &cik) in maxwell[i].iter().enumerate().take(n) {
            row_sum += cik;
            if i != k {
                mutual[i][k] = -cik;
            }
        }
        mutual[i][i] = row_sum;
    }
    mutual
}

/// Result of dual-box `1/R` open-boundary Maxwell capacitance matrix extrapolation (`Mode G2`, `FR-31`),
/// with Dykstra projection onto the nearest Stieltjes cone (symmetric Z-matrix with `C_ij <= 0` for `i != j`,
/// non-negative row sums `\sum_j C_ij >= 0`, and strictly SPD) when `0 < raw_offdiag_pos_ratio <= 0.05`.
#[derive(Clone, Debug, PartialEq)]
pub struct OpenBoxExtrapolation {
    pub raw_maxwell_f: Vec<Vec<f64>>,
    pub projected_stieltjes_f: Vec<Vec<f64>>,
    pub mutual_matrix_f: Vec<Vec<f64>>,
    pub raw_min_eig: f64,
    pub projected_min_eig: f64,
    pub raw_offdiag_pos_ratio: f64,
    pub projected_offdiag_pos_ratio: f64,
    pub used_projection: bool,
    pub note: String,
}

impl OpenBoxExtrapolation {
    pub fn to_json(&self) -> Value {
        json!({
            "raw_maxwell_F": self.raw_maxwell_f,
            "projected_stieltjes_F": self.projected_stieltjes_f,
            "mutual_matrix_F": self.mutual_matrix_f,
            "raw_min_eig": self.raw_min_eig,
            "projected_min_eig": self.projected_min_eig,
            "raw_offdiag_pos_ratio": self.raw_offdiag_pos_ratio,
            "projected_offdiag_pos_ratio": self.projected_offdiag_pos_ratio,
            "used_projection": self.used_projection,
            "note": self.note,
        })
    }
}

fn maxwell_min_eig_and_offdiag_pos_ratio(mat: &[Vec<f64>]) -> Result<(f64, f64), String> {
    let n = mat.len();
    if n == 0 {
        return Err("capacitance matrix must be non-empty".into());
    }
    let mut flat = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            flat.push(0.5 * (mat[i][j] + mat[j][i]));
        }
    }
    let (eigs, _) = dense_eigh(&flat, n).map_err(|e| format!("dense_eigh failed: {e}"))?;
    let min_eig = eigs.iter().copied().fold(f64::INFINITY, f64::min);
    let mut max_pos_ratio = 0.0_f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let denom = (mat[i][i].abs() * mat[j][j].abs()).sqrt().max(1.0e-30);
            let ratio = (0.5 * (mat[i][j] + mat[j][i])) / denom;
            if ratio > max_pos_ratio {
                max_pos_ratio = ratio;
            }
        }
    }
    Ok((min_eig, max_pos_ratio))
}

fn project_onto_stieltjes_cone(raw: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = raw.len();
    let max_diag = (0..n).map(|i| raw[i][i].abs()).fold(0.0_f64, f64::max).max(1.0e-30);
    let eig_floor = 1.0e-12 * max_diag;
    let diag_floor = 1.0e-12 * max_diag;

    let mut x = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            x[i][j] = 0.5 * (raw[i][j] + raw[j][i]);
        }
    }

    // Dykstra's cyclic alternating projection onto:
    //   0: Z-matrix cone { X = X^T : X_ij <= 0 (i != j), X_ii >= diag_floor }
    //   1..=n: Row-sum half-spaces H_i = { X = X^T : \sum_j X_ij >= 0 }
    //   n+1: Spectral SPD cone { X = X^T : X >= eig_floor * I }
    let mut p_z = vec![vec![0.0; n]; n];
    let mut p_rows = vec![vec![vec![0.0; n]; n]; n];
    let mut p_spd = vec![vec![0.0; n]; n];
    let inv_norm_hi_sq = 2.0 / ((n + 1) as f64);

    for _iter in 0..40 {
        // 1. Project onto Z-matrix cone
        let mut y = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                y[i][j] = x[i][j] + p_z[i][j];
            }
        }
        let mut x_next = vec![vec![0.0; n]; n];
        for i in 0..n {
            x_next[i][i] = y[i][i].max(diag_floor);
            for j in (i + 1)..n {
                let sym = 0.5 * (y[i][j] + y[j][i]);
                let clamped = sym.min(0.0);
                x_next[i][j] = clamped;
                x_next[j][i] = clamped;
            }
        }
        for i in 0..n {
            for j in 0..n {
                p_z[i][j] = y[i][j] - x_next[i][j];
            }
        }
        x = x_next;

        // 2. Project onto each symmetric row-sum half-space H_i (\sum_j X_ij >= 0)
        for r in 0..n {
            let mut yr = vec![vec![0.0; n]; n];
            for i in 0..n {
                for j in 0..n {
                    yr[i][j] = x[i][j] + p_rows[r][i][j];
                }
            }
            let row_sum: f64 = (0..n).map(|j| 0.5 * (yr[r][j] + yr[j][r])).sum();
            let mut xr = yr.clone();
            if row_sum < 0.0 {
                let alpha = -row_sum * inv_norm_hi_sq;
                xr[r][r] += alpha;
                for j in 0..n {
                    if j != r {
                        xr[r][j] += 0.5 * alpha;
                        xr[j][r] += 0.5 * alpha;
                    }
                }
            }
            for i in 0..n {
                for j in 0..n {
                    p_rows[r][i][j] = yr[i][j] - xr[i][j];
                }
            }
            x = xr;
        }

        // 3. Project onto spectral SPD cone via Higham eigendecomposition
        let mut ys = vec![vec![0.0; n]; n];
        let mut flat = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                let v = 0.5 * ((x[i][j] + p_spd[i][j]) + (x[j][i] + p_spd[j][i]));
                ys[i][j] = v;
                flat.push(v);
            }
        }
        let (eigs, vecs) = dense_eigh(&flat, n).map_err(|e| format!("dense_eigh failed: {e}"))?;
        let mut xs = vec![vec![0.0; n]; n];
        for (k, &lam) in eigs.iter().enumerate() {
            let lk = lam.max(eig_floor);
            let vk = &vecs[k];
            for i in 0..n {
                for j in i..n {
                    let prod = lk * vk[i] * vk[j];
                    xs[i][j] += prod;
                    if i != j {
                        xs[j][i] += prod;
                    }
                }
            }
        }
        for i in 0..n {
            for j in 0..n {
                p_spd[i][j] = ys[i][j] - xs[i][j];
            }
        }
        x = xs;
    }

    // Final exact boundary polish guaranteeing X_ij <= 0 (i != j) and \sum_j X_ij >= 0 (hence diagonally dominant SPD)
    for i in 0..n {
        for j in (i + 1)..n {
            let sym = (0.5 * (x[i][j] + x[j][i])).min(0.0);
            x[i][j] = sym;
            x[j][i] = sym;
        }
    }
    for i in 0..n {
        let off_sum: f64 = (0..n).filter(|&j| j != i).map(|j| -x[i][j]).sum();
        if x[i][i] < off_sum + eig_floor {
            x[i][i] = off_sum + eig_floor;
        }
    }
    Ok(x)
}

/// Computes dual-box `1/R` open-boundary Maxwell capacitance extrapolation (`Mode G2`, `FR-31`):
/// `C_inf = (R_2 * C_2 - R_1 * C_1) / (R_2 - R_1)` symmetrized, evaluates `raw_offdiag_pos_ratio`,
/// and applies Dykstra projection onto the nearest Stieltjes cone when `0 < raw_offdiag_pos_ratio <= 0.05`.
pub fn extrapolate_open_box_maxwell(
    c1: &[Vec<f64>],
    r1: f64,
    c2: &[Vec<f64>],
    r2: f64,
) -> Result<OpenBoxExtrapolation, String> {
    let n = c1.len();
    if n == 0 {
        return Err("extrapolate_open_box_maxwell requires non-empty capacitance matrices".into());
    }
    if c2.len() != n || c1.iter().any(|r| r.len() != n) || c2.iter().any(|r| r.len() != n) {
        return Err(format!("c1 and c2 must both be {n}x{n} square matrices"));
    }
    if c1.iter().flatten().any(|x| !x.is_finite()) || c2.iter().flatten().any(|x| !x.is_finite()) {
        return Err("c1 and c2 must contain only finite capacitance values".into());
    }
    if !(r1.is_finite() && r2.is_finite() && r1 > 0.0 && r2 > 0.0) {
        return Err(format!("box radii r1 ({r1:e}) and r2 ({r2:e}) must be finite and positive"));
    }
    let dr = r2 - r1;
    if dr.abs() <= 1.0e-15 * r1.max(r2) {
        return Err(format!("box radii r1 ({r1:e}) and r2 ({r2:e}) must be distinct"));
    }

    let mut raw = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i..n {
            let v_ij = (r2 * c2[i][j] - r1 * c1[i][j]) / dr;
            let v_ji = (r2 * c2[j][i] - r1 * c1[j][i]) / dr;
            let sym = 0.5 * (v_ij + v_ji);
            raw[i][j] = sym;
            raw[j][i] = sym;
        }
    }

    let (raw_min_eig, raw_offdiag_pos_ratio) = maxwell_min_eig_and_offdiag_pos_ratio(&raw)?;
    let has_neg_row_sum = (0..n).any(|i| raw[i].iter().sum::<f64>() < -1.0e-15 * raw[i][i].abs());
    let needs_stieltjes_proj = (raw_offdiag_pos_ratio > 0.0 && raw_offdiag_pos_ratio <= 0.05)
        || (raw_offdiag_pos_ratio <= 0.05 && (raw_min_eig <= 0.0 || has_neg_row_sum));

    let (projected_stieltjes_f, used_projection) =
        if needs_stieltjes_proj { (project_onto_stieltjes_cone(&raw)?, true) } else { (raw.clone(), false) };
    let (projected_min_eig, projected_offdiag_pos_ratio) =
        maxwell_min_eig_and_offdiag_pos_ratio(&projected_stieltjes_f)?;
    let mutual_matrix_f = mutual_from_maxwell(&projected_stieltjes_f);

    let note = if used_projection {
        format!(
            "Applied Dykstra Stieltjes cone projection to remove dual-box 1/R discretization overshoot (raw_offdiag_pos_ratio = {:.3e} -> {:.3e}); prefer single-box far_field(order = 2) (AF-3) for single-solve M-matrix preservation.",
            raw_offdiag_pos_ratio, projected_offdiag_pos_ratio
        )
    } else if raw_offdiag_pos_ratio > 0.05 {
        format!(
            "Raw dual-box 1/R extrapolation exceeded Stieltjes projection threshold (raw_offdiag_pos_ratio = {:.3e} > 0.05); prefer single-box far_field(order = 2) (AF-3).",
            raw_offdiag_pos_ratio
        )
    } else {
        "Raw dual-box 1/R extrapolation already satisfies Stieltjes M-matrix invariants; prefer single-box far_field(order = 2) (AF-3) to avoid dual-mesh solve overhead.".to_string()
    };

    Ok(OpenBoxExtrapolation {
        raw_maxwell_f: raw,
        projected_stieltjes_f,
        mutual_matrix_f,
        raw_min_eig,
        projected_min_eig,
        raw_offdiag_pos_ratio,
        projected_offdiag_pos_ratio,
        used_projection,
        note,
    })
}

/// Eigendecomposes symmetric PSD `maxwell` (`m x m`), returning `(ec_matrix_hz, condition_number, null_modes)`.
#[allow(clippy::needless_range_loop)]
fn spectral_ec_and_cond(maxwell: &[Vec<f64>]) -> Result<(Vec<Vec<f64>>, f64, usize), Error> {
    let m = maxwell.len();
    if m == 0 {
        return Ok((Vec::new(), 1.0, 0));
    }
    let mut flat = Vec::with_capacity(m * m);
    for row in maxwell {
        flat.extend_from_slice(row);
    }
    let (eigs, vecs) = dense_eigh(&flat, m).map_err(Error::from_kernel)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    let lam_min = eigs.iter().copied().fold(f64::INFINITY, f64::min);
    if lam_max <= 0.0 || lam_min < -1.0e-8 * lam_max {
        return Err(Error::invalid("reduced Maxwell capacitance matrix is not positive semi-definite"));
    }

    let null_thresh = NULL_MODE_REL_TOL * lam_max;
    let null_modes = eigs.iter().filter(|&&lam| lam <= null_thresh).count();
    let denom = lam_min.max(lam_max / COND_CAP);
    let condition_number = (lam_max / denom).min(COND_CAP);

    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    let mut ec = vec![vec![0.0; m]; m];
    for (k, &lam) in eigs.iter().enumerate().take(m) {
        if lam > null_thresh {
            let inv_l = scale / lam;
            let vk = &vecs[k];
            for i in 0..m {
                for j in 0..m {
                    ec[i][j] += inv_l * vk[i] * vk[j];
                }
            }
        }
    }
    for i in 0..m {
        for j in (i + 1)..m {
            let sym = 0.5 * (ec[i][j] + ec[j][i]);
            ec[i][j] = sym;
            ec[j][i] = sym;
        }
    }
    Ok((ec, condition_number, null_modes))
}

/// Asserts that `maxwell` (`n x n`) is strictly symmetric positive definite (SPD), returning its `ec_matrix_Hz`.
fn assert_spd_and_ec(maxwell: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, Error> {
    let n = maxwell.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut flat = Vec::with_capacity(n * n);
    for row in maxwell {
        flat.extend_from_slice(row);
    }
    let (eigs, vecs) = dense_eigh(&flat, n).map_err(Error::from_kernel)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    let lam_min = eigs.iter().copied().fold(f64::INFINITY, f64::min);
    if !(lam_max > 0.0 && lam_min > 1.0e-12 * lam_max) {
        return Err(Error::invalid(format!(
            "thickness-corrected Maxwell capacitance matrix is not strictly SPD (lambda_min={lam_min:e}, lambda_max={lam_max:e})"
        )));
    }
    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    let mut ec = vec![vec![0.0; n]; n];
    for (k, &lam) in eigs.iter().enumerate() {
        let inv_l = scale / lam;
        let vk = &vecs[k];
        for i in 0..n {
            for j in 0..n {
                ec[i][j] += inv_l * vk[i] * vk[j];
            }
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let sym = 0.5 * (ec[i][j] + ec[j][i]);
            ec[i][j] = sym;
            ec[j][i] = sym;
        }
    }
    Ok(ec)
}

/// Default calibrated 2D Schwarz–Christoffel matched-asymptotic coefficients `(c0, c1)` for a planar
/// rectangular step of thickness `t > 0` in medium `eps_vac` atop substrate `eps_sub`.
///
/// In homogeneous space (`eps_sub == eps_vac`), exact Schwarz–Christoffel mapping of a semi-infinite plate
/// of thickness `t` combined with Eshelby's crack-tip J-integral gives `c1 = 1 / (2 pi)` and
/// `c0 = (1 + ln(4 pi * alpha_s)) / (2 pi)`. On a dielectric substrate (`eps_sub > eps_vac`) where the
/// metal film protrudes into the upper `eps_vac` half-space while `\tilde{I}_{edge}` integrates across
/// both `eps_vac` and `eps_sub`, the outer singularity intensity scales as `(eps_vac + eps_sub) / 2` while
/// the sidewall capacitance enhancement in `eps_vac` refracts into the substrate with interface enhancement
/// `eta(eps_sub / eps_vac)`.
fn default_sc_thickness_coeffs(eps_vac: f64, eps_sub: f64) -> (f64, f64) {
    let ev = eps_vac.max(1.0);
    let es = eps_sub.max(ev);
    // Exact dielectric-step scaling relative to `(ev + es) / 2` outer edge-tube normalization:
    // For a metal film of thickness `t` protruding into the upper `ev` half-space (`z in [0, t]`) above a flat
    // substrate (`z < 0`, permittivity `es`), `\tilde{I}_{edge}` scales with `(ev + es) / 2` while the vertical
    // step sidewall sits in `ev` enhanced by positive substrate image refraction `1 + 0.146 * (es - ev)/(es + ev)`.
    let contrast = (es - ev) / (es + ev);
    let r_dielec = (2.0 * ev / (ev + es)) * (1.0 + 0.146 * contrast);
    let c1 = (1.0 / (2.0 * PI)) * r_dielec;
    let c0 = 0.450 * r_dielec;
    (c0, c1)
}

pub(crate) struct ThicknessCorrectionOut {
    pub(crate) maxwell_f: Vec<Vec<f64>>,
    pub(crate) mutual_f: Vec<Vec<f64>>,
    pub(crate) ec_hz: Vec<Vec<f64>>,
    pub(crate) info_json: Value,
}

#[derive(Clone, Debug)]
pub(crate) struct PreparedThicknessGram {
    pub(crate) i_edge: Vec<Vec<f64>>,
    pub(crate) s_eff: f64,
    pub(crate) eps_vac: f64,
    pub(crate) eps_sub: f64,
    pub(crate) r0: f64,
    pub(crate) def_c0: f64,
    pub(crate) def_c1: f64,
    pub(crate) c0_override: Option<f64>,
    pub(crate) c1_override: Option<f64>,
    pub(crate) perimeter_segments: usize,
    pub(crate) perimeter_len_m: f64,
}

/// Integrates the geometry-only 3D matched-asymptotic edge-tube Gram matrix `\tilde{I}_{edge}` once so zero-solve
/// UQ can vary `metal_thickness_m`, `c0`, and `c1` without re-integrating the 3D mesh.
pub(crate) fn prepare_3d_thickness_gram(
    mesh: &Mesh,
    result: &SolveResult,
    cap: &CapacitanceMatrixData,
    params_map: &Map<String, Value>,
    t_metal_nom: f64,
) -> Result<PreparedThicknessGram, Error> {
    if !(t_metal_nom.is_finite() && t_metal_nom > 0.0) {
        return Err(Error::invalid("parameter 'metal_thickness_m' must be a positive finite thickness in m"));
    }
    if result.dim != 3 {
        return Err(Error::invalid(format!(
            "metal_thickness_m requires a 3D volumetric electrostatic result, got {}D",
            result.dim
        )));
    }
    if result.log.iter().any(|line| line.contains("planar_bem")) {
        return Err(Error::invalid(
            "metal_thickness_m requires 3D volumetric_fem electric fields; planar_bem stores boundary-only dummy fields",
        ));
    }

    let edge_cutoff_opt = match params_map.get("edge_cutoff_m") {
        None | Some(Value::Null) => None,
        Some(v) => Some(
            v.as_f64()
                .filter(|&x| x.is_finite() && x > 0.0)
                .ok_or_else(|| Error::invalid("parameter 'edge_cutoff_m' must be a positive finite radius in m"))?,
        ),
    };
    let c0_override = match params_map.get("c0") {
        None | Some(Value::Null) => None,
        Some(v) => Some(
            v.as_f64()
                .filter(|x| x.is_finite())
                .ok_or_else(|| Error::invalid("parameter 'c0' must be a finite number"))?,
        ),
    };
    let c1_override = match params_map.get("c1") {
        None | Some(Value::Null) => None,
        Some(v) => Some(
            v.as_f64()
                .filter(|x| x.is_finite())
                .ok_or_else(|| Error::invalid("parameter 'c1' must be a finite number"))?,
        ),
    };

    let disc = discretise(mesh)?;
    if disc.dim != 3 {
        return Err(Error::invalid(format!("metal_thickness_m requires a 3D mesh, got {}D", disc.dim)));
    }
    let nc = disc.n_cells();
    let n_term = cap.terminals.len();
    if result.solutions.len() < n_term {
        return Err(Error::invalid(format!(
            "metal_thickness_m requires {} terminal electric-field solution(s) (found {}; ensure save_fields=true and volumetric_fem)",
            n_term,
            result.solutions.len()
        )));
    }

    let mut accessors = Vec::with_capacity(n_term);
    for i in 0..n_term {
        accessors.push(EFieldAccessor::new(&result.solutions[i], 3, nc)?);
    }
    let order = accessors[0].order;

    let meta = TlsMeta::of(result).unwrap_or_default();
    let cell_eps: Vec<f64> =
        (0..nc).map(|ci| meta.domain_eps_r.get(&disc.cell_labels[ci]).copied().unwrap_or(1.0)).collect();
    let geo = geometry(&disc);
    let quad3 = quad_rule(3, 4);
    let shapes3: Vec<Vec<f64>> = quad3.lam.iter().map(|l| shape(l, order).0).collect();

    // Guard against dummy boundary-only fields (e.g. stripped-log planar_bem):
    // In planar_bem, V is 1 on terminal boundary vertices and 0 on all free interior vertices, so every cell
    // disjoint from terminal vertices has identically zero electric field AND volumetric energy != C_{00}.
    let mut c00_vol = 0.0;
    let mut non_zero_cells = 0usize;
    for ci in 0..nc {
        let mut cell_e2 = 0.0;
        for (q, sh) in shapes3.iter().enumerate() {
            let e = accessors[0].eval_real_with_shapes(ci, sh);
            cell_e2 += quad3.w[q] * (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        }
        if cell_e2 > 0.0 {
            non_zero_cells += 1;
        }
        c00_vol += EPS0 * cell_eps[ci] * geo.vol[ci] * cell_e2;
    }
    if c00_vol <= 0.0 || (nc > 16 && non_zero_cells * 4 < nc) {
        return Err(Error::invalid(
            "metal_thickness_m rejected dummy/incomplete electric fields (expected full 3D volumetric_fem solution)",
        ));
    }

    let perim = detect_sheet_perimeter_3d(&disc, &meta, &cell_eps);
    if perim.bvh.is_empty() {
        return Err(Error::invalid(
            "metal_thickness_m found no interior zero-thickness sheet perimeter edges on the 3D mesh",
        ));
    }

    let r0 = edge_cutoff_opt.unwrap_or_else(|| (1.15 * perim.h_edge_mean).max(4.0 * t_metal_nom));
    let r0_sq = r0 * r0;
    let sub_tet = subdivided_tet_quad(4);
    let sub_tet_shapes: Vec<Vec<f64>> = sub_tet.iter().map(|(l, _)| shape(l, order).0).collect();

    // Integrate n_term x n_term edge-tube Gram matrix:
    // I_tilde_{ij} = (1 / R_0) * int_{r <= R_0} eps_0 * eps_r(x) * E^(i)(x) . E^(j)(x) dV
    let mut i_edge = vec![vec![0.0; n_term]; n_term];
    let mut e_at_q = vec![[0.0; 3]; n_term];

    for ci in 0..nc {
        let v = disc.cell(ci);
        let pts = [pt3(&disc, v[0]), pt3(&disc, v[1]), pt3(&disc, v[2]), pt3(&disc, v[3])];
        let mut d_min = f64::INFINITY;
        let mut d_max = 0.0_f64;
        for &pt in &pts {
            let d = perim.bvh.min_dist(pt);
            d_min = d_min.min(d);
            d_max = d_max.max(d);
        }
        let h_tet = (6.0 * geo.vol[ci]).cbrt();
        if d_min > r0 + h_tet {
            continue;
        }
        let vol_scale = (EPS0 * cell_eps[ci] * geo.vol[ci]) / r0;
        if d_max <= r0 {
            for (q, sh) in shapes3.iter().enumerate() {
                let w = vol_scale * quad3.w[q];
                for t_idx in 0..n_term {
                    e_at_q[t_idx] = accessors[t_idx].eval_real_with_shapes(ci, sh);
                }
                for ia in 0..n_term {
                    let ea = e_at_q[ia];
                    for ib in ia..n_term {
                        let eb = e_at_q[ib];
                        i_edge[ia][ib] += w * (ea[0] * eb[0] + ea[1] * eb[1] + ea[2] * eb[2]);
                    }
                }
            }
        } else {
            for (q, &(l4, w_ref)) in sub_tet.iter().enumerate() {
                let x = [
                    l4[0] * pts[0][0] + l4[1] * pts[1][0] + l4[2] * pts[2][0] + l4[3] * pts[3][0],
                    l4[0] * pts[0][1] + l4[1] * pts[1][1] + l4[2] * pts[2][1] + l4[3] * pts[3][1],
                    l4[0] * pts[0][2] + l4[1] * pts[1][2] + l4[2] * pts[2][2] + l4[3] * pts[3][2],
                ];
                if perim.bvh.min_dist_sq(x) <= r0_sq {
                    let w = vol_scale * w_ref;
                    let sh = &sub_tet_shapes[q];
                    for t_idx in 0..n_term {
                        e_at_q[t_idx] = accessors[t_idx].eval_real_with_shapes(ci, sh);
                    }
                    for ia in 0..n_term {
                        let ea = e_at_q[ia];
                        for ib in ia..n_term {
                            let eb = e_at_q[ib];
                            i_edge[ia][ib] += w * (ea[0] * eb[0] + ea[1] * eb[1] + ea[2] * eb[2]);
                        }
                    }
                }
            }
        }
    }
    for ia in 0..n_term {
        for ib in (ia + 1)..n_term {
            i_edge[ib][ia] = i_edge[ia][ib];
        }
    }

    let tr_cm: f64 = (0..n_term).map(|i| cap.maxwell_matrix_f[i][i]).sum();
    let tr_iedge: f64 = (0..n_term).map(|i| i_edge[i][i]).sum();
    if !(tr_iedge > 0.0 && tr_cm > 0.0) {
        return Err(Error::invalid("edge-tube Gram matrix trace is zero; cannot evaluate s_eff"));
    }
    let s_eff = tr_cm / tr_iedge;
    if s_eff <= t_metal_nom {
        return Err(Error::invalid(format!(
            "metal_thickness_m ({t_metal_nom:e} m) exceeds effective gap scale s_eff ({s_eff:e} m); asymptotic thin-film regime violated"
        )));
    }

    let eps_vac = cell_eps.iter().copied().fold(f64::INFINITY, f64::min).max(1.0);
    let eps_sub = cell_eps.iter().copied().fold(eps_vac, f64::max);
    let (def_c0, def_c1) = default_sc_thickness_coeffs(eps_vac, eps_sub);

    Ok(PreparedThicknessGram {
        i_edge,
        s_eff,
        eps_vac,
        eps_sub,
        r0,
        def_c0,
        def_c1,
        c0_override,
        c1_override,
        perimeter_segments: perim.bvh.len(),
        perimeter_len_m: perim.total_perim_len,
    })
}

/// Evaluates the 3D thickness capacitance correction from a cached [`PreparedThicknessGram`].
pub(crate) fn apply_3d_thickness_gram(
    cap: &CapacitanceMatrixData,
    gram: &PreparedThicknessGram,
    t_metal: f64,
    c0_override: Option<f64>,
    c1_override: Option<f64>,
) -> Result<ThicknessCorrectionOut, Error> {
    if !(t_metal.is_finite() && t_metal > 0.0) {
        return Err(Error::invalid("parameter 'metal_thickness_m' must be a positive finite thickness in m"));
    }
    if gram.s_eff <= t_metal {
        return Err(Error::invalid(format!(
            "metal_thickness_m ({t_metal:e} m) exceeds effective gap scale s_eff ({:e} m); asymptotic thin-film regime violated",
            gram.s_eff
        )));
    }
    let n_term = cap.terminals.len();
    let c0 = c0_override.or(gram.c0_override).unwrap_or(gram.def_c0);
    let c1 = c1_override.or(gram.c1_override).unwrap_or(gram.def_c1);

    let factor = t_metal * (c0 + c1 * (gram.s_eff / t_metal).ln());
    let mut delta_c = vec![vec![0.0; n_term]; n_term];
    let mut maxwell_corr = vec![vec![0.0; n_term]; n_term];
    for ia in 0..n_term {
        for ib in ia..n_term {
            let base_sym = 0.5 * (cap.maxwell_matrix_f[ia][ib] + cap.maxwell_matrix_f[ib][ia]);
            let d_ij = factor * gram.i_edge[ia][ib];
            delta_c[ia][ib] = d_ij;
            delta_c[ib][ia] = d_ij;
            maxwell_corr[ia][ib] = base_sym + d_ij;
            maxwell_corr[ib][ia] = base_sym + d_ij;
        }
    }

    let ec_corr = assert_spd_and_ec(&maxwell_corr)?;
    let mutual_corr = mutual_from_maxwell(&maxwell_corr);

    let info_json = json!({
        "metal_thickness_m": t_metal,
        "edge_cutoff_m": gram.r0,
        "s_eff_m": gram.s_eff,
        "c0": c0,
        "c1": c1,
        "eps_vac": gram.eps_vac,
        "eps_sub": gram.eps_sub,
        "perimeter_segments": gram.perimeter_segments,
        "perimeter_len_m": gram.perimeter_len_m,
        "i_edge_F_per_m": gram.i_edge,
        "delta_maxwell_F": delta_c,
        "maxwell_uncorrected_F": cap.maxwell_matrix_f,
    });

    Ok(ThicknessCorrectionOut { maxwell_f: maxwell_corr, mutual_f: mutual_corr, ec_hz: ec_corr, info_json })
}

/// Computes the 3D matched-asymptotic edge-tube Gram-matrix capacitance thickness correction (FR-11 Pt 1).
fn compute_3d_thickness_correction(
    mesh: &Mesh,
    result: &SolveResult,
    cap: &CapacitanceMatrixData,
    params_map: &Map<String, Value>,
    t_metal: f64,
) -> Result<ThicknessCorrectionOut, Error> {
    let gram = prepare_3d_thickness_gram(mesh, result, cap, params_map, t_metal)?;
    apply_3d_thickness_gram(cap, &gram, t_metal, None, None)
}

fn parse_ground_names(params_map: &Map<String, Value>) -> Result<Vec<String>, Error> {
    match params_map.get("ground") {
        None | Some(Value::Null) => Ok(Vec::new()),
        Some(Value::String(s)) => Ok(vec![s.clone()]),
        Some(Value::Array(arr)) => {
            let mut out = Vec::with_capacity(arr.len());
            for (i, item) in arr.iter().enumerate() {
                let s = item
                    .as_str()
                    .ok_or_else(|| Error::invalid(format!("ground[{i}] must be a terminal name string")))?;
                out.push(s.to_string());
            }
            Ok(out)
        }
        _ => Err(Error::invalid("parameter 'ground' must be a terminal name string or array of strings")),
    }
}

fn parse_eliminate_mode(params_map: &Map<String, Value>) -> Result<EliminateMode, Error> {
    match params_map.get("eliminate") {
        None | Some(Value::Null) => Ok(EliminateMode::Drop),
        Some(Value::String(s)) => match s.as_str() {
            "drop" => Ok(EliminateMode::Drop),
            "floating_schur" => Ok(EliminateMode::FloatingSchur),
            other => {
                Err(Error::invalid(format!("parameter 'eliminate' must be 'drop' or 'floating_schur', not '{other}'")))
            }
        },
        _ => Err(Error::invalid("parameter 'eliminate' must be 'drop' or 'floating_schur'")),
    }
}

fn reduce_maxwell_matrix(
    raw: &[Vec<f64>],
    sig_indices: &[usize],
    ground_indices: &[usize],
    eliminate: EliminateMode,
    strict_psd_gg: bool,
) -> Result<Vec<Vec<f64>>, Error> {
    let m = sig_indices.len();
    let mut reduced = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in a..m {
            let (ia, ib) = (sig_indices[a], sig_indices[b]);
            let sym = 0.5 * (raw[ia][ib] + raw[ib][ia]);
            reduced[a][b] = sym;
            reduced[b][a] = sym;
        }
    }

    if eliminate == EliminateMode::FloatingSchur && !ground_indices.is_empty() {
        let kg = ground_indices.len();
        let mut flat_gg = Vec::with_capacity(kg * kg);
        for &ig in ground_indices {
            for &jg in ground_indices {
                flat_gg.push(0.5 * (raw[ig][jg] + raw[jg][ig]));
            }
        }
        let (eigs_gg, vecs_gg) = dense_eigh(&flat_gg, kg).map_err(Error::from_kernel)?;
        let lam_max_gg = eigs_gg.iter().copied().fold(0.0_f64, f64::max);
        let lam_min_gg = eigs_gg.iter().copied().fold(f64::INFINITY, f64::min);
        if (strict_psd_gg && lam_max_gg <= 0.0) || lam_min_gg < -1.0e-8 * lam_max_gg.max(1.0e-30) {
            return Err(Error::invalid("ground submatrix C_gg is not positive semi-definite in floating_schur"));
        }
        if lam_max_gg > 0.0 {
            let thresh_gg = NULL_MODE_REL_TOL * lam_max_gg;
            let mut c_sg = vec![vec![0.0; kg]; m];
            for (a, &ia) in sig_indices.iter().enumerate() {
                for (p, &ig) in ground_indices.iter().enumerate() {
                    c_sg[a][p] = 0.5 * (raw[ia][ig] + raw[ig][ia]);
                }
            }
            for (k, &lam) in eigs_gg.iter().enumerate().take(kg) {
                if lam > thresh_gg {
                    let vk = &vecs_gg[k];
                    let w: Vec<f64> = (0..m).map(|a| (0..kg).map(|p| c_sg[a][p] * vk[p]).sum()).collect();
                    for a in 0..m {
                        for b in a..m {
                            let sub = (w[a] * w[b]) / lam;
                            reduced[a][b] -= sub;
                            if a != b {
                                reduced[b][a] -= sub;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(reduced)
}

fn qubit_c_sigma_from_psd_bound(c_red: &[Vec<f64>], pad_indices: &[usize]) -> Option<f64> {
    if pad_indices.len() == 1 {
        let c_ss = c_red[pad_indices[0]][pad_indices[0]];
        return (c_ss.is_finite() && c_ss >= 0.0).then_some(c_ss);
    }
    let (s, b) = (pad_indices[0], pad_indices[1]);
    let c_ss = c_red[s][s];
    let c_bb = c_red[b][b];
    let c_sb = 0.5 * (c_red[s][b] + c_red[b][s]);
    if !(c_ss.is_finite() && c_bb.is_finite() && c_sb.is_finite()) || c_ss < 0.0 || c_bb < 0.0 {
        return None;
    }
    let sum_diag = c_ss + c_bb;
    if sum_diag <= 0.0 {
        return Some(0.0);
    }
    let denom = sum_diag + 2.0 * c_sb;
    let det = (c_ss * c_bb - c_sb * c_sb).max(0.0);
    if denom > 1.0e-14 * sum_diag {
        Some(det / denom)
    } else if denom >= -1.0e-12 * sum_diag && c_sb <= 0.0 {
        Some((-c_sb).max(0.0))
    } else {
        None
    }
}

type MatrixBoundsSlice<'a> = (&'a [Vec<f64>], &'a [Vec<f64>]);

#[allow(clippy::too_many_arguments)]
fn parse_qubits(
    q_val: &Value,
    all_terminals: &[String],
    ground_set: &BTreeSet<String>,
    reduced_terminals: &[String],
    red_idx_of: &BTreeMap<&str, usize>,
    c_red: &[Vec<f64>],
    c_red_bounds: Option<MatrixBoundsSlice<'_>>,
) -> Result<Vec<ComputedQubit>, Error> {
    let arr = q_val.as_array().ok_or_else(|| Error::invalid("parameter 'qubits' must be an array of objects"))?;
    let mut seen_names = BTreeSet::new();
    let mut out = Vec::with_capacity(arr.len());
    let ec_scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);

    for (qi, item) in arr.iter().enumerate() {
        let obj = item.as_object().ok_or_else(|| Error::invalid(format!("qubits[{qi}] must be an object")))?;
        for k in obj.keys() {
            if !matches!(k.as_str(), "name" | "pads" | "E_J_Hz") {
                return Err(Error::invalid(format!(
                    "unknown field '{k}' in qubits[{qi}] (expected: name, pads, E_J_Hz)"
                )));
            }
        }
        let name = obj
            .get("name")
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| Error::invalid(format!("qubits[{qi}].name must be a non-empty string")))?
            .to_string();
        if !seen_names.insert(name.clone()) {
            return Err(Error::invalid(format!("duplicate qubit name '{name}'")));
        }

        let pads_arr = obj
            .get("pads")
            .and_then(Value::as_array)
            .ok_or_else(|| Error::invalid(format!("qubit '{name}' requires 'pads' as [str] or [str, str]")))?;
        if !(1..=2).contains(&pads_arr.len()) {
            return Err(Error::invalid(format!(
                "qubit '{name}' 'pads' must have 1 (grounded) or 2 (differential) terminal names, got {}",
                pads_arr.len()
            )));
        }
        let mut pads = Vec::with_capacity(pads_arr.len());
        let mut pad_indices = Vec::with_capacity(pads_arr.len());
        for (pi, pv) in pads_arr.iter().enumerate() {
            let p_name =
                pv.as_str().ok_or_else(|| Error::invalid(format!("qubit '{name}' pads[{pi}] must be a string")))?;
            if !all_terminals.iter().any(|t| t == p_name) {
                return Err(Error::invalid(format!(
                    "qubit '{name}' pad '{p_name}' does not match any terminal in {all_terminals:?}"
                )));
            }
            if ground_set.contains(p_name) {
                return Err(Error::invalid(format!(
                    "qubit '{name}' pad '{p_name}' cannot also be an eliminated ground terminal"
                )));
            }
            let &idx = red_idx_of.get(p_name).ok_or_else(|| {
                Error::invalid(format!("qubit '{name}' pad '{p_name}' is not present in reduced terminals"))
            })?;
            pads.push(p_name.to_string());
            pad_indices.push(idx);
        }
        if pad_indices.len() == 2 && pad_indices[0] == pad_indices[1] {
            return Err(Error::invalid(format!("qubit '{name}' two-pad specification must name two distinct pads")));
        }

        let e_j_hz = match obj.get("E_J_Hz") {
            None | Some(Value::Null) => None,
            Some(v) => Some(v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                Error::invalid(format!("qubit '{name}' E_J_Hz must be a finite positive frequency in Hz"))
            })?),
        };

        let mut c_eff_map = BTreeMap::new();
        let mut beta_map = BTreeMap::new();
        let (topology, c_sigma_f, extra_pad_caps) = if pad_indices.len() == 1 {
            let s = pad_indices[0];
            let c_ss = c_red[s][s];
            if !(c_ss > 0.0 && c_ss.is_finite()) {
                return Err(Error::invalid(format!("qubit '{name}' self-capacitance C_ss must be positive")));
            }
            for (k_idx, term_k) in reduced_terminals.iter().enumerate() {
                if k_idx == s {
                    continue;
                }
                let c_eff_k = c_red[s][k_idx].abs();
                c_eff_map.insert(term_k.clone(), c_eff_k);
                beta_map.insert(term_k.clone(), c_eff_k / c_ss);
            }
            ("grounded_1pad", c_ss, None)
        } else {
            let (s, b) = (pad_indices[0], pad_indices[1]);
            let c_ss = c_red[s][s];
            let c_bb = c_red[b][b];
            let c_sb = c_red[s][b];
            if !(c_ss > 0.0 && c_bb > 0.0) {
                return Err(Error::invalid(format!("qubit '{name}' pad self-capacitances must be positive")));
            }
            let sum_diag = c_ss + c_bb;
            let denom = sum_diag + 2.0 * c_sb;
            let det = c_ss * c_bb - c_sb * c_sb;
            let c_sigma = if denom > 1.0e-14 * sum_diag {
                if det <= 0.0 {
                    return Err(Error::invalid(format!(
                        "qubit '{name}' 2-pad Maxwell submatrix is not positive definite"
                    )));
                }
                det / denom
            } else if denom >= -1.0e-12 * sum_diag && c_sb < 0.0 {
                -c_sb
            } else {
                return Err(Error::invalid(format!(
                    "qubit '{name}' 2-pad Maxwell entries (C_ss={c_ss:e}, C_bb={c_bb:e}, C_sb={c_sb:e}) violate passivity"
                )));
            };

            for (k_idx, term_k) in reduced_terminals.iter().enumerate() {
                if k_idx == s || k_idx == b {
                    continue;
                }
                let c_sk = c_red[s][k_idx].abs();
                let c_bk = c_red[b][k_idx].abs();
                let c_eff_k = if denom > 1.0e-14 * sum_diag {
                    ((c_bb + c_sb) * c_sk - (c_ss + c_sb) * c_bk).abs() / denom
                } else {
                    0.5 * (c_sk - c_bk).abs()
                };
                c_eff_map.insert(term_k.clone(), c_eff_k);
                beta_map.insert(term_k.clone(), c_eff_k / c_sigma);
            }
            ("differential_2pad", c_sigma, Some((c_ss, c_bb, c_sb)))
        };

        let e_c_hz = ec_scale / c_sigma_f;
        let mut q_json = Map::new();
        q_json.insert("name".into(), json!(name));
        q_json.insert("topology".into(), json!(topology));
        q_json.insert("pads".into(), json!(pads));
        if let Some((c_ss, c_bb, c_sb)) = extra_pad_caps {
            q_json.insert("C_ss_F".into(), json!(c_ss));
            q_json.insert("C_bb_F".into(), json!(c_bb));
            q_json.insert("C_sb_F".into(), json!(c_sb));
        }
        q_json.insert("C_Sigma_F".into(), json!(c_sigma_f));
        q_json.insert("E_C_Hz".into(), json!(e_c_hz));
        if let Some((c_low_red, c_up_red)) = c_red_bounds {
            if let (Some(cs_low_raw), Some(cs_up_raw)) = (
                qubit_c_sigma_from_psd_bound(c_low_red, &pad_indices),
                qubit_c_sigma_from_psd_bound(c_up_red, &pad_indices),
            ) {
                if cs_up_raw > 0.0 {
                    if cs_low_raw > cs_up_raw * (1.0 + 1.0e-8) {
                        return Err(Error::invalid(format!(
                            "qubit '{name}' lower bound C_Sigma_low_F ({cs_low_raw:e}) exceeds upper bound C_Sigma_up_F ({cs_up_raw:e})"
                        )));
                    }
                    let cs_low = cs_low_raw.min(cs_up_raw);
                    let cs_up = cs_up_raw.max(cs_low_raw);
                    let ec_low = ec_scale / cs_up;
                    q_json.insert("C_Sigma_low_F".into(), json!(cs_low));
                    q_json.insert("C_Sigma_up_F".into(), json!(cs_up));
                    q_json.insert("E_C_low_Hz".into(), json!(ec_low));
                    if cs_low > 0.0 {
                        q_json.insert("E_C_up_Hz".into(), json!(ec_scale / cs_low));
                    }
                }
            }
        }
        q_json.insert("C_eff_F".into(), json!(c_eff_map));
        q_json.insert("beta".into(), json!(beta_map));

        let (mut f01_opt, mut nzpf_opt) = (None, None);
        if let Some(ej) = e_j_hz {
            let ej_over_ec = ej / e_c_hz;
            let f01 = (8.0 * ej * e_c_hz).sqrt() - e_c_hz;
            let alpha = -e_c_hz;
            let nzpf = (ej / (32.0 * e_c_hz)).powf(0.25);
            q_json.insert("E_J_Hz".into(), json!(ej));
            q_json.insert("ej_over_ec".into(), json!(ej_over_ec));
            q_json.insert("f01_Hz".into(), json!(f01));
            q_json.insert("anharmonicity_Hz".into(), json!(alpha));
            q_json.insert("n_zpf".into(), json!(nzpf));
            f01_opt = Some(f01);
            nzpf_opt = Some(nzpf);
        }

        out.push(ComputedQubit {
            name,
            pads,
            c_sigma_f,
            e_c_hz,
            c_eff_f: c_eff_map,
            beta: beta_map,
            e_j_hz,
            f01_hz: f01_opt,
            n_zpf: nzpf_opt,
            json: Value::Object(q_json),
        });
    }
    Ok(out)
}

fn parse_resonators(
    r_val: &Value,
    all_terminals: &[String],
    ground_set: &BTreeSet<String>,
    red_idx_of: &BTreeMap<&str, usize>,
) -> Result<Vec<ComputedResonator>, Error> {
    let arr = r_val.as_array().ok_or_else(|| Error::invalid("parameter 'resonators' must be an array of objects"))?;
    let mut seen_names = BTreeSet::new();
    let mut out = Vec::with_capacity(arr.len());
    let hbar = PLANCK_H / (2.0 * PI);

    for (ri, item) in arr.iter().enumerate() {
        let obj = item.as_object().ok_or_else(|| Error::invalid(format!("resonators[{ri}] must be an object")))?;
        for k in obj.keys() {
            if !matches!(k.as_str(), "name" | "terminal" | "f_r_Hz" | "Z0_ohm" | "resonator_type" | "C_r_F") {
                return Err(Error::invalid(format!(
                    "unknown field '{k}' in resonators[{ri}] (expected: name, terminal, f_r_Hz, Z0_ohm, resonator_type, C_r_F)"
                )));
            }
        }
        let name = obj
            .get("name")
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| Error::invalid(format!("resonators[{ri}].name must be a non-empty string")))?
            .to_string();
        if !seen_names.insert(name.clone()) {
            return Err(Error::invalid(format!("duplicate resonator name '{name}'")));
        }

        let terminal = obj
            .get("terminal")
            .and_then(Value::as_str)
            .ok_or_else(|| Error::invalid(format!("resonator '{name}' requires 'terminal' string")))?
            .to_string();
        if !all_terminals.iter().any(|t| t == &terminal) {
            return Err(Error::invalid(format!(
                "resonator '{name}' terminal '{terminal}' does not match any terminal in {all_terminals:?}"
            )));
        }
        if ground_set.contains(&terminal) || !red_idx_of.contains_key(terminal.as_str()) {
            return Err(Error::invalid(format!(
                "resonator '{name}' terminal '{terminal}' cannot be an eliminated ground terminal"
            )));
        }

        let f_r_hz = obj
            .get("f_r_Hz")
            .and_then(Value::as_f64)
            .filter(|&x| x.is_finite() && x > 0.0)
            .ok_or_else(|| Error::invalid(format!("resonator '{name}' requires finite positive 'f_r_Hz'")))?;

        let z0_ohm = match obj.get("Z0_ohm") {
            None | Some(Value::Null) => 50.0,
            Some(v) => v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                Error::invalid(format!("resonator '{name}' 'Z0_ohm' must be a finite positive impedance in Ohm"))
            })?,
        };

        let resonator_type = match obj.get("resonator_type") {
            None | Some(Value::Null) => "quarter_wave",
            Some(Value::String(s)) => match s.as_str() {
                "quarter_wave" | "half_wave" | "lumped" => s.as_str(),
                other => {
                    return Err(Error::invalid(format!(
                        "resonator '{name}' 'resonator_type' must be 'quarter_wave', 'half_wave', or 'lumped', not '{other}'"
                    )))
                }
            },
            _ => return Err(Error::invalid(format!("resonator '{name}' 'resonator_type' must be a string"))),
        };

        let c_r_override = match obj.get("C_r_F") {
            None | Some(Value::Null) => None,
            Some(v) => Some(v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                Error::invalid(format!("resonator '{name}' 'C_r_F' must be a finite positive capacitance in F"))
            })?),
        };

        let omega_r = 2.0 * PI * f_r_hz;
        let c_r_f = match resonator_type {
            "quarter_wave" => PI / (4.0 * omega_r * z0_ohm),
            "half_wave" => PI / (2.0 * omega_r * z0_ohm),
            "lumped" => c_r_override.ok_or_else(|| {
                Error::invalid(format!("resonator '{name}' with resonator_type='lumped' requires positive 'C_r_F'"))
            })?,
            _ => unreachable!(),
        };
        let v_rms_v = (hbar * omega_r / (2.0 * c_r_f)).sqrt();

        let r_json = json!({
            "name": name,
            "terminal": terminal,
            "resonator_type": resonator_type,
            "f_r_Hz": f_r_hz,
            "Z0_ohm": z0_ohm,
            "C_r_F": c_r_f,
            "V_rms_V": v_rms_v,
        });
        out.push(ComputedResonator { name, terminal, f_r_hz, v_rms_v, json: r_json });
    }
    Ok(out)
}

fn insert_optional_bound_matrices(out: &mut Map<String, Value>, cap: &CapacitanceMatrixData) {
    if let Some(low) = &cap.maxwell_lower_matrix_f {
        out.insert("maxwell_lower_matrix_F".into(), json!(low));
    }
    if let Some(up) = &cap.maxwell_upper_matrix_f {
        out.insert("maxwell_upper_matrix_F".into(), json!(up));
    }
    if let Some(mid) = &cap.maxwell_mid_matrix_f {
        out.insert("maxwell_mid_matrix_F".into(), json!(mid));
    }
    if let Some(err) = &cap.maxwell_error_bound_f {
        out.insert("maxwell_error_bound_F".into(), json!(err));
    }
    if let Some(cond) = &cap.conductance_matrix_s {
        out.insert("conductance_matrix_S".into(), json!(cond));
    }
    if let Some(qd) = &cap.q_dielectric_matrix {
        out.insert("q_dielectric_matrix".into(), json!(qd));
    }
    if !cap.adapt_history.is_empty() {
        out.insert("adapt_history".into(), json!(cap.adapt_history));
    }
}

fn validate_optional_bound_matrix(name: &str, mat_opt: &Option<Vec<Vec<f64>>>, n: usize) -> Result<(), Error> {
    if let Some(mat) = mat_opt {
        if mat.len() != n || mat.iter().any(|r| r.len() != n) {
            return Err(Error::invalid(format!("{name} dimensions do not match terminal count {n}")));
        }
        if mat.iter().flatten().any(|x| !x.is_finite()) {
            return Err(Error::NonFinite(format!("{name} contains non-finite entries")));
        }
    }
    Ok(())
}

fn run_capacitance_matrix_impl(mesh: Option<&Mesh>, result: &SolveResult, params: &Value) -> Result<Value, Error> {
    let cap = capacitance_matrix(result)?;
    let base_reply = |cap: &CapacitanceMatrixData| -> Value {
        let mut out = Map::new();
        out.insert("terminals".into(), json!(cap.terminals));
        out.insert("maxwell_matrix_F".into(), json!(cap.maxwell_matrix_f));
        out.insert("mutual_matrix_F".into(), json!(cap.mutual_matrix_f));
        out.insert("ec_matrix_Hz".into(), json!(cap.ec_matrix_hz));
        insert_optional_bound_matrices(&mut out, cap);
        Value::Object(out)
    };
    let params_map = match params {
        Value::Null => return Ok(base_reply(&cap)),
        Value::Object(m) if m.is_empty() => return Ok(base_reply(&cap)),
        Value::Object(m) => m,
        _ => return Err(Error::invalid("capacitance_matrix params must be a JSON object or null")),
    };

    for k in params_map.keys() {
        if !matches!(
            k.as_str(),
            "ground" | "eliminate" | "qubits" | "resonators" | "metal_thickness_m" | "edge_cutoff_m" | "c0" | "c1"
        ) {
            return Err(Error::invalid(format!(
                "unknown capacitance_matrix parameter '{k}' (expected: ground, eliminate, qubits, resonators, metal_thickness_m, edge_cutoff_m, c0, c1)"
            )));
        }
    }

    let n = cap.terminals.len();
    if n == 0 || cap.maxwell_matrix_f.len() != n || cap.maxwell_matrix_f.iter().any(|r| r.len() != n) {
        return Err(Error::invalid("CapacitanceMatrixData dimensions are empty or inconsistent"));
    }
    if cap.maxwell_matrix_f.iter().flatten().any(|x| !x.is_finite()) {
        return Err(Error::NonFinite("Maxwell capacitance matrix contains non-finite entries".into()));
    }
    validate_optional_bound_matrix("maxwell_lower_matrix_f", &cap.maxwell_lower_matrix_f, n)?;
    validate_optional_bound_matrix("maxwell_upper_matrix_f", &cap.maxwell_upper_matrix_f, n)?;
    validate_optional_bound_matrix("maxwell_mid_matrix_f", &cap.maxwell_mid_matrix_f, n)?;
    validate_optional_bound_matrix("maxwell_error_bound_f", &cap.maxwell_error_bound_f, n)?;

    // Optional FR-11 Pt 1 3D matched-asymptotic finite-thickness correction BEFORE Schur/qubit reduction
    let thickness_corr = match params_map.get("metal_thickness_m") {
        None | Some(Value::Null) => {
            if params_map.get("edge_cutoff_m").is_some_and(|v| !v.is_null())
                || params_map.get("c0").is_some_and(|v| !v.is_null())
                || params_map.get("c1").is_some_and(|v| !v.is_null())
            {
                return Err(Error::invalid(
                    "parameters 'edge_cutoff_m', 'c0', and 'c1' require 'metal_thickness_m' to be specified",
                ));
            }
            None
        }
        Some(v) => {
            let t = v
                .as_f64()
                .ok_or_else(|| Error::invalid("parameter 'metal_thickness_m' must be a positive finite number (m)"))?;
            let m_ref = mesh.ok_or_else(|| {
                Error::invalid(
                    "metal_thickness_m requires the 3D Mesh; call capacitance_matrix_with_mesh_and_params or CapacitanceMatrixAnalysis::run_with_params",
                )
            })?;
            Some(compute_3d_thickness_correction(m_ref, result, &cap, params_map, t)?)
        }
    };

    run_capacitance_matrix_from_thickness(&cap, thickness_corr, params_map)
}

/// Evaluates Schur/ground reduction, qubit extraction, and resonator couplings from `cap` and an optional
/// precomputed 3D thickness correction (used directly by zero-solve UQ without re-integrating the 3D mesh).
pub(crate) fn run_capacitance_matrix_from_thickness(
    cap: &CapacitanceMatrixData,
    thickness_corr: Option<ThicknessCorrectionOut>,
    params_map: &Map<String, Value>,
) -> Result<Value, Error> {
    let active_maxwell = thickness_corr.as_ref().map_or(&cap.maxwell_matrix_f, |tc| &tc.maxwell_f);
    let active_mutual = thickness_corr.as_ref().map_or(&cap.mutual_matrix_f, |tc| &tc.mutual_f);
    let active_ec = thickness_corr.as_ref().map_or(&cap.ec_matrix_hz, |tc| &tc.ec_hz);

    let n = cap.terminals.len();
    let ground_raw = parse_ground_names(params_map)?;
    let eliminate = parse_eliminate_mode(params_map)?;
    let mut ground_idx_set = BTreeSet::new();
    for g in &ground_raw {
        let idx = cap.terminals.iter().position(|t| t == g).ok_or_else(|| {
            Error::invalid(format!("ground terminal '{g}' does not match any terminal in {:?}", cap.terminals))
        })?;
        ground_idx_set.insert(idx);
    }

    let ground_indices: Vec<usize> = ground_idx_set.iter().copied().collect();
    let ground_terminals: Vec<String> = ground_indices.iter().map(|&i| cap.terminals[i].clone()).collect();
    let ground_name_set: BTreeSet<String> = ground_terminals.iter().cloned().collect();
    let sig_indices: Vec<usize> = (0..n).filter(|i| !ground_idx_set.contains(i)).collect();
    if sig_indices.is_empty() {
        return Err(Error::invalid("ground cannot eliminate all terminals; at least one signal terminal must remain"));
    }

    let reduced_terminals: Vec<String> = sig_indices.iter().map(|&i| cap.terminals[i].clone()).collect();
    let maxwell_reduced_f = reduce_maxwell_matrix(active_maxwell, &sig_indices, &ground_indices, eliminate, true)?;
    let reduced_bounds = match (thickness_corr.is_none(), &cap.maxwell_lower_matrix_f, &cap.maxwell_upper_matrix_f) {
        (true, Some(low), Some(up)) => {
            let red_low = reduce_maxwell_matrix(low, &sig_indices, &ground_indices, eliminate, false)?;
            let red_up = reduce_maxwell_matrix(up, &sig_indices, &ground_indices, eliminate, true)?;
            Some((red_low, red_up))
        }
        _ => None,
    };

    let mutual_reduced_f = mutual_from_maxwell(&maxwell_reduced_f);
    let (ec_reduced_hz, condition_number, null_modes) = spectral_ec_and_cond(&maxwell_reduced_f)?;
    let has_null_mode = null_modes > 0;

    let red_idx_of: BTreeMap<&str, usize> =
        reduced_terminals.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();

    let qubits_opt = match params_map.get("qubits") {
        None | Some(Value::Null) => None,
        Some(v) => Some(parse_qubits(
            v,
            &cap.terminals,
            &ground_name_set,
            &reduced_terminals,
            &red_idx_of,
            &maxwell_reduced_f,
            reduced_bounds.as_ref().map(|(lo, hi)| (lo.as_slice(), hi.as_slice())),
        )?),
    };

    let resonators_opt = match params_map.get("resonators") {
        None | Some(Value::Null) => None,
        Some(v) => Some(parse_resonators(v, &cap.terminals, &ground_name_set, &red_idx_of)?),
    };

    let mut out = Map::new();
    out.insert("terminals".into(), json!(cap.terminals));
    out.insert("maxwell_matrix_F".into(), json!(active_maxwell));
    out.insert("mutual_matrix_F".into(), json!(active_mutual));
    out.insert("ec_matrix_Hz".into(), json!(active_ec));
    if thickness_corr.is_none() {
        insert_optional_bound_matrices(&mut out, cap);
    }
    out.insert("eliminate".into(), json!(eliminate.as_str()));
    out.insert("ground_terminals".into(), json!(ground_terminals));
    out.insert("reduced_terminals".into(), json!(reduced_terminals));
    out.insert("maxwell_reduced_F".into(), json!(maxwell_reduced_f));
    if let Some((red_low, red_up)) = &reduced_bounds {
        out.insert("maxwell_reduced_lower_F".into(), json!(red_low));
        out.insert("maxwell_reduced_upper_F".into(), json!(red_up));
    }
    out.insert("mutual_reduced_F".into(), json!(mutual_reduced_f));
    out.insert("ec_reduced_Hz".into(), json!(ec_reduced_hz));
    out.insert("condition_number".into(), json!(condition_number));
    out.insert("null_modes".into(), json!(null_modes));
    out.insert("has_null_mode".into(), json!(has_null_mode));

    if let Some(tc) = thickness_corr {
        out.insert("thickness_correction".into(), tc.info_json);
    }

    if let (Some(qubits), Some(resonators)) = (&qubits_opt, &resonators_opt) {
        let mut couplings = Vec::with_capacity(qubits.len() * resonators.len());
        for q in qubits {
            for r in resonators {
                if q.pads.iter().any(|p| p == &r.terminal) {
                    return Err(Error::invalid(format!(
                        "resonator '{}' terminal '{}' cannot be one of qubit '{}' pads {:?}",
                        r.name, r.terminal, q.name, q.pads
                    )));
                }
                let c_eff = q.c_eff_f.get(&r.terminal).copied().unwrap_or(0.0);
                let beta = q.beta.get(&r.terminal).copied().unwrap_or_else(|| {
                    if q.c_sigma_f > 0.0 {
                        c_eff / q.c_sigma_f
                    } else {
                        0.0
                    }
                });
                let mut c_obj = Map::new();
                c_obj.insert("qubit".into(), json!(q.name));
                c_obj.insert("resonator".into(), json!(r.name));
                c_obj.insert("terminal".into(), json!(r.terminal));
                c_obj.insert("C_eff_F".into(), json!(c_eff));
                c_obj.insert("beta".into(), json!(beta));
                c_obj.insert("V_rms_V".into(), json!(r.v_rms_v));
                if let (Some(_ej), Some(f01), Some(nzpf)) = (q.e_j_hz, q.f01_hz, q.n_zpf) {
                    let g_hz = 2.0 * E_CHARGE * beta * nzpf * r.v_rms_v / PLANCK_H;
                    let delta_hz = f01 - r.f_r_hz;
                    let chi_hz = if delta_hz.abs() > 1.0e-12 * r.f_r_hz.max(1.0)
                        && (delta_hz - q.e_c_hz).abs() > 1.0e-12 * q.e_c_hz.max(1.0)
                    {
                        Some(-(g_hz * g_hz) / delta_hz * (q.e_c_hz / (delta_hz - q.e_c_hz)))
                    } else {
                        None
                    };
                    c_obj.insert("n_zpf".into(), json!(nzpf));
                    c_obj.insert("g_Hz".into(), json!(g_hz));
                    c_obj.insert("Delta_Hz".into(), json!(delta_hz));
                    c_obj.insert("chi_Hz".into(), json!(chi_hz));
                }
                couplings.push(Value::Object(c_obj));
            }
        }
        out.insert("couplings".into(), Value::Array(couplings));
    }

    if let Some(qubits) = qubits_opt {
        out.insert("qubits".into(), Value::Array(qubits.into_iter().map(|q| q.json).collect()));
    }
    if let Some(resonators) = resonators_opt {
        out.insert("resonators".into(), Value::Array(resonators.into_iter().map(|r| r.json).collect()));
    }

    Ok(Value::Object(out))
}

/// Runs the multi-terminal capacitance analysis without a mesh (`params` for ground elimination, qubit
/// reduction, and resonator dispersive coupling).
///
/// When `params` is `Value::Null` or an empty object `{}`, returns the exact 4-field base JSON
/// (`terminals`, `maxwell_matrix_F`, `mutual_matrix_F`, `ec_matrix_Hz`) unchanged.
pub fn capacitance_matrix_with_params(result: &SolveResult, params: &Value) -> Result<Value, Error> {
    run_capacitance_matrix_impl(None, result, params)
}

/// Runs the multi-terminal capacitance analysis with the 3D [`Mesh`] (enabling `metal_thickness_m`
/// matched-asymptotic `I_edge` Gram-matrix thickness correction prior to Schur/qubit reduction).
pub fn capacitance_matrix_with_mesh_and_params(
    mesh: &Mesh,
    result: &SolveResult,
    params: &Value,
) -> Result<Value, Error> {
    run_capacitance_matrix_impl(Some(mesh), result, params)
}

/// Multi-terminal Maxwell, mutual, and charging-energy capacitance matrix analysis behind the [`Analysis`] trait.
pub struct CapacitanceMatrixAnalysis;

impl Analysis for CapacitanceMatrixAnalysis {
    fn id(&self) -> &'static str {
        "capacitance_matrix"
    }

    fn name(&self) -> &'static str {
        "Multi-terminal capacitance and charging energy matrix"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        capacitance_matrix_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        capacitance_matrix_with_params(result, &Value::Null)
    }

    fn run_with_params(&self, mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        capacitance_matrix_with_mesh_and_params(mesh, result, params)
    }
}

/// Extracts the typed [`InductanceMatrixData`] payload from `result`.
pub fn inductance_matrix(result: &SolveResult) -> Result<InductanceMatrixData, Error> {
    InductanceMatrixData::of(result)
        .ok_or_else(|| Error::invalid("the result carries no inductance matrix (run a magnetostatic study of dfs)"))
}

/// True when `result` carries a magnetostatic inductance matrix payload.
pub fn inductance_matrix_applies(result: &SolveResult) -> bool {
    InductanceMatrixData::of(result).is_some()
}

/// Multi-terminal magnetostatic inductance matrix analysis behind the [`Analysis`] trait.
pub struct InductanceMatrixAnalysis;

impl Analysis for InductanceMatrixAnalysis {
    fn id(&self) -> &'static str {
        "inductance_matrix"
    }

    fn name(&self) -> &'static str {
        "Multi-terminal inductance matrix"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        inductance_matrix_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        let ind = inductance_matrix(result)?;
        Ok(json!({
            "terminals": ind.terminals,
            "inductance_matrix_H": ind.inductance_matrix_h,
            "formulation": ind.formulation,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capacitance_and_inductance_analyses_round_trip_to_json() {
        let sm = femcore::testkit::rect(1.0, 0.5, 2, 2);
        let mesh = Mesh::from_solver_mesh(&sm).unwrap();

        let cap = CapacitanceMatrixData {
            terminals: vec!["pad1".into(), "pad2".into()],
            maxwell_matrix_f: vec![vec![80e-15, -30e-15], vec![-30e-15, 85e-15]],
            mutual_matrix_f: vec![vec![50e-15, 30e-15], vec![30e-15, 55e-15]],
            ec_matrix_hz: vec![vec![2.8e8, 1.0e8], vec![1.0e8, 2.6e8]],
            ..Default::default()
        };
        let ind = InductanceMatrixData {
            terminals: vec!["loop1".into()],
            inductance_matrix_h: vec![vec![150e-12]],
            formulation: "mixed_3d".into(),
        };
        let res = SolveResult {
            physics: "dfs".into(),
            study: "electrostatic".into(),
            scalars: vec![cap.scalar(), ind.scalar()],
            ..Default::default()
        };

        assert!(CapacitanceMatrixAnalysis.applies_to(&res));
        assert!(InductanceMatrixAnalysis.applies_to(&res));

        let cap_json = CapacitanceMatrixAnalysis.run(&mesh, &res).unwrap();
        assert_eq!(cap_json["terminals"], json!(["pad1", "pad2"]));
        assert_eq!(cap_json["maxwell_matrix_F"][0][0].as_f64().unwrap(), 80e-15);
        assert_eq!(cap_json["mutual_matrix_F"][0][1].as_f64().unwrap(), 30e-15);
        assert_eq!(cap_json["ec_matrix_Hz"][1][1].as_f64().unwrap(), 2.6e8);

        let ind_json = InductanceMatrixAnalysis.run(&mesh, &res).unwrap();
        assert_eq!(ind_json["terminals"], json!(["loop1"]));
        assert_eq!(ind_json["inductance_matrix_H"][0][0].as_f64().unwrap(), 150e-12);
        assert_eq!(ind_json["formulation"], json!("mixed_3d"));

        let empty = SolveResult::default();
        assert!(!CapacitanceMatrixAnalysis.applies_to(&empty));
        assert!(!InductanceMatrixAnalysis.applies_to(&empty));
        assert!(CapacitanceMatrixAnalysis.run(&mesh, &empty).is_err());
        assert!(InductanceMatrixAnalysis.run(&mesh, &empty).is_err());
    }

    #[test]
    fn fr31_extrapolate_open_box_maxwell_projects_small_overshoot_onto_stieltjes_cone() {
        let r1 = 550.0e-6;
        let r2 = 850.0e-6;

        // Case 1: Clean dual-box 1/R matrices already satisfying Stieltjes M-matrix invariants
        let c1_clean = vec![
            vec![95.0e-15, -6.2e-15, -72.0e-15],
            vec![-6.2e-15, 48.0e-15, -34.0e-15],
            vec![-72.0e-15, -34.0e-15, 290.0e-15],
        ];
        let c2_clean = vec![
            vec![92.0e-15, -6.1e-15, -70.0e-15],
            vec![-6.1e-15, 46.5e-15, -33.2e-15],
            vec![-70.0e-15, -33.2e-15, 240.0e-15],
        ];
        let clean = extrapolate_open_box_maxwell(&c1_clean, r1, &c2_clean, r2).unwrap();
        assert!(!clean.used_projection);
        assert_eq!(clean.raw_offdiag_pos_ratio, 0.0);
        assert_eq!(clean.projected_offdiag_pos_ratio, 0.0);
        assert_eq!(clean.raw_maxwell_f, clean.projected_stieltjes_f);
        assert!(clean.projected_min_eig > 0.0);
        assert!(clean.note.contains("far_field(order = 2)"));

        // Case 2: Weakly coupled pair (0, 1) where C_01(R1) = -1.20 fF and C_01(R2) = -0.58 fF
        // so raw 1/R extrapolation C_01,inf = (850 * (-0.58) - 550 * (-1.20)) / 300 = +0.5567 fF > 0
        // (`raw_offdiag_pos_ratio ~ 8.3e-3` in (0, 0.05], modeling s1_1_p2_00_r0148 Mode G2 overshoot)
        let mut c1_overshoot = c1_clean.clone();
        let mut c2_overshoot = c2_clean.clone();
        c1_overshoot[0][1] = -1.20e-15;
        c1_overshoot[1][0] = -1.20e-15;
        c2_overshoot[0][1] = -0.58e-15;
        c2_overshoot[1][0] = -0.58e-15;

        let proj = extrapolate_open_box_maxwell(&c1_overshoot, r1, &c2_overshoot, r2).unwrap();
        assert!(proj.used_projection, "expected Stieltjes projection for 0 < raw_offdiag_pos_ratio <= 0.05");
        assert!(proj.raw_offdiag_pos_ratio > 1.0e-3 && proj.raw_offdiag_pos_ratio <= 0.05);
        assert_eq!(proj.projected_offdiag_pos_ratio, 0.0);
        assert!(proj.projected_min_eig > 0.0);
        for i in 0..3 {
            let row_sum: f64 = proj.projected_stieltjes_f[i].iter().sum();
            assert!(row_sum >= 0.0, "row {i} sum must be non-negative: {row_sum:e}");
            for j in 0..3 {
                assert!(
                    (proj.projected_stieltjes_f[i][j] - proj.projected_stieltjes_f[j][i]).abs() <= 1.0e-28,
                    "projected matrix must be symmetric"
                );
                if i != j {
                    assert!(
                        proj.projected_stieltjes_f[i][j] <= 0.0,
                        "off-diagonal C[{i},{j}] must be <= 0, got {:e}",
                        proj.projected_stieltjes_f[i][j]
                    );
                }
            }
        }
        let j_val = proj.to_json();
        assert_eq!(j_val["used_projection"], json!(true));
        assert!(j_val["note"].as_str().unwrap().contains("far_field(order = 2)"));

        // Case 3: Validation errors on invalid radii or dimensions
        assert!(extrapolate_open_box_maxwell(&c1_clean, r1, &c2_clean, r1).is_err());
        assert!(extrapolate_open_box_maxwell(&c1_clean, -r1, &c2_clean, r2).is_err());
        assert!(extrapolate_open_box_maxwell(&[], r1, &[], r2).is_err());
    }
}
