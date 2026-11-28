// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Multi-port Hermite Tangential Galerkin projection-based reduced-order modeling (PROM)
//! with isometric QR image-space residual certification (`ImageProm`), shared between
//! electromagnetic (`nyx-ems`) and piezoelectric (`nyx-sms`) driven frequency sweeps.
//!
//! At each greedy anchor $\omega_a$, a single sparse factorization solves for the port fields
//! $\mathbf{X}(\omega_a)$ and their first two frequency derivatives $\mathbf{X}'(\omega_a)$,
//! $\mathbf{X}''(\omega_a)$. After optional block-partitioned SVD truncation ($s_k \ge 10^{-6} s_0$),
//! orthonormal real basis vectors build the Galerkin subspace $\mathbf{Q}$ while $\hat{\mathbf{A}}_m \mathbf{Q}$
//! and $\hat{\mathbf{b}}_q$ build an isometric real image basis $\mathbf{U}$ for microsecond
//! full-space residual evaluation across all sweep frequencies.
#![allow(clippy::needless_range_loop)]

use crate::linalg::blr::{BlrOptions, ShiftedBlrLu};
use crate::linalg::csr::{c, CCoo, CCsr, C};
use crate::linalg::direct::{FactorBackend, ShiftedSolver};
use crate::linalg::driver::{backward_residual_many, thin_svd_c64, DirectSolveCheck, DirectSolver};
use std::ops::Range;
use std::time::Instant;

/// `[re, im]` dense port matrix table (`P x P`).
pub type PortMatrixPairs = Vec<Vec<[f64; 2]>>;

/// Euclidean inner product of two equal-length real slices.
#[inline]
pub fn dot_real(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

/// Non-conjugate dot product $\mathbf{a}^T \mathbf{b}$ of a real slice and a complex slice.
#[inline]
pub fn dot_real_complex(a: &[f64], b: &[C]) -> C {
    let mut acc = c(0.0, 0.0);
    for (&x, &y) in a.iter().zip(b) {
        acc += y * x;
    }
    acc
}

/// Maps `FactorBackend` to a concise canonical label matching driven solver logs.
pub fn format_factor_backend(b: FactorBackend) -> &'static str {
    match b {
        FactorBackend::None => "None",
        FactorBackend::CholeskyAmd => "CholeskyAmd",
        FactorBackend::CholeskyNd => "CholeskyNd",
        FactorBackend::LbltAmd => "LbltAmd",
        FactorBackend::LbltNd => "LbltNd",
        FactorBackend::LuColamd => "LuColamd",
        FactorBackend::LuAmd => "LuAmd",
        FactorBackend::LuNd => "LuNd",
        FactorBackend::BlrNd => "BlrNd",
    }
}

/// Two-pass Modified Gram-Schmidt $l_2$-orthonormalization of `cand` into `basis`.
pub fn append_l2_orthonormal(basis: &mut Vec<Vec<f64>>, cand: &[f64]) -> bool {
    let norm0 = dot_real(cand, cand).sqrt();
    if norm0 <= 1e-30 || !norm0.is_finite() {
        return false;
    }
    let mut v = cand.to_vec();
    for _pass in 0..2 {
        for q in basis.iter() {
            let proj = dot_real(&v, q);
            for (vi, &qi) in v.iter_mut().zip(q) {
                *vi -= proj * qi;
            }
        }
    }
    let norm = dot_real(&v, &v).sqrt();
    if norm <= 1e-12 * norm0 || !norm.is_finite() {
        return false;
    }
    let inv = 1.0 / norm;
    for vi in &mut v {
        *vi *= inv;
    }
    basis.push(v);
    true
}

/// Solves an $r \times r$ complex dense linear system $A X = B$ for $P$ right-hand sides (`b_cols[q]` of length $r$)
/// using Gaussian elimination with partial pivoting.
pub fn solve_dense_complex_multi(a_mat: &[C], b_cols: &[Vec<C>], r: usize) -> Option<Vec<Vec<C>>> {
    let p = b_cols.len();
    if r == 0 {
        return Some(vec![Vec::new(); p]);
    }
    let mut a = a_mat.to_vec();
    let mut rhs = vec![c(0.0, 0.0); r * p];
    for q in 0..p {
        for i in 0..r {
            rhs[i * p + q] = b_cols[q][i];
        }
    }

    for k in 0..r {
        let mut piv_row = k;
        let mut piv_max = a[k * r + k].norm();
        for i in (k + 1)..r {
            let val = a[i * r + k].norm();
            if val > piv_max {
                piv_max = val;
                piv_row = i;
            }
        }
        if piv_max <= 1e-300 || !piv_max.is_finite() {
            return None;
        }
        if piv_row != k {
            for j in k..r {
                a.swap(k * r + j, piv_row * r + j);
            }
            for q in 0..p {
                rhs.swap(k * p + q, piv_row * p + q);
            }
        }
        let inv_piv = c(1.0, 0.0) / a[k * r + k];
        for i in (k + 1)..r {
            let factor = a[i * r + k] * inv_piv;
            if factor.re == 0.0 && factor.im == 0.0 {
                continue;
            }
            for j in (k + 1)..r {
                let akj = a[k * r + j];
                a[i * r + j] -= factor * akj;
            }
            for q in 0..p {
                let rkq = rhs[k * p + q];
                rhs[i * p + q] -= factor * rkq;
            }
        }
    }

    let mut x_cols = vec![vec![c(0.0, 0.0); r]; p];
    for i in (0..r).rev() {
        let inv_diag = c(1.0, 0.0) / a[i * r + i];
        for q in 0..p {
            let mut sum = rhs[i * p + q];
            for j in (i + 1)..r {
                sum -= a[i * r + j] * x_cols[q][j];
            }
            x_cols[q][i] = sum * inv_diag;
        }
    }
    Some(x_cols)
}

/// Converts a $P \times P$ scattering matrix $\mathbf{S}$ and reference impedances $Z_{0,p}$ into
/// `(S, Z, Y)` arrays stored as `[re, im]` pairs, using
/// $\mathbf{Z} = \mathbf{Z}_0^{1/2}(\mathbf{I}+\mathbf{S})(\mathbf{I}-\mathbf{S})^{-1}\mathbf{Z}_0^{1/2}$
/// and $\mathbf{Y} = \mathbf{Z}_0^{-1/2}(\mathbf{I}-\mathbf{S})(\mathbf{I}+\mathbf{S})^{-1}\mathbf{Z}_0^{-1/2}$.
pub fn s_to_z_and_y(
    s_mat: &[Vec<C>],
    z0_ports: &[f64],
) -> (PortMatrixPairs, PortMatrixPairs, PortMatrixPairs) {
    let (s, z, y, _rep) = s_to_z_and_y_checked(s_mat, z0_ports);
    (s, z, y)
}

/// Diagnostic metrics produced alongside $S \leftrightarrow (Z, Y)$ port conversions and
/// scattering-matrix physicality checks.
#[derive(Debug, Clone, Default)]
pub struct NetworkConversionReport {
    pub used_regularization: bool,
    pub singular_fallback: bool,
    pub max_solve_residual: f64,
    pub max_pivot_ratio: f64,
    pub warning_messages: Vec<String>,
}

/// Solves $A X = B$ for $P$ right-hand sides and returns `(x_cols, rel_backward_residual, pivot_ratio)`.
fn solve_dense_complex_multi_checked(
    a_mat: &[C],
    b_cols: &[Vec<C>],
    r: usize,
) -> Option<(Vec<Vec<C>>, f64, f64)> {
    let p = b_cols.len();
    if r == 0 {
        return Some((vec![Vec::new(); p], 0.0, 1.0));
    }
    let mut a = a_mat.to_vec();
    let mut rhs = vec![c(0.0, 0.0); r * p];
    for q in 0..p {
        for i in 0..r {
            rhs[i * p + q] = b_cols[q][i];
        }
    }

    let mut min_piv = f64::INFINITY;
    let mut max_piv = 0.0_f64;
    for k in 0..r {
        let mut piv_row = k;
        let mut piv_max = a[k * r + k].norm();
        for i in (k + 1)..r {
            let val = a[i * r + k].norm();
            if val > piv_max {
                piv_max = val;
                piv_row = i;
            }
        }
        min_piv = min_piv.min(piv_max);
        max_piv = max_piv.max(piv_max);
        if piv_max <= 1e-300 || !piv_max.is_finite() {
            return None;
        }
        if piv_row != k {
            for j in k..r {
                a.swap(k * r + j, piv_row * r + j);
            }
            for q in 0..p {
                rhs.swap(k * p + q, piv_row * p + q);
            }
        }
        let inv_piv = c(1.0, 0.0) / a[k * r + k];
        for i in (k + 1)..r {
            let factor = a[i * r + k] * inv_piv;
            if factor.re == 0.0 && factor.im == 0.0 {
                continue;
            }
            for j in (k + 1)..r {
                let akj = a[k * r + j];
                a[i * r + j] -= factor * akj;
            }
            for q in 0..p {
                let rkq = rhs[k * p + q];
                rhs[i * p + q] -= factor * rkq;
            }
        }
    }

    let mut x_cols = vec![vec![c(0.0, 0.0); r]; p];
    for i in (0..r).rev() {
        let inv_diag = c(1.0, 0.0) / a[i * r + i];
        for q in 0..p {
            let mut sum = rhs[i * p + q];
            for j in (i + 1)..r {
                sum -= a[i * r + j] * x_cols[q][j];
            }
            x_cols[q][i] = sum * inv_diag;
        }
    }

    let a_fro = a_mat.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
    let mut worst_res = 0.0_f64;
    for q in 0..p {
        let x_norm = x_cols[q].iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
        let b_norm = b_cols[q].iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
        let mut r2 = 0.0_f64;
        for i in 0..r {
            let mut ax_i = c(0.0, 0.0);
            for j in 0..r {
                ax_i += a_mat[i * r + j] * x_cols[q][j];
            }
            r2 += (ax_i - b_cols[q][i]).norm_sqr();
        }
        let res = r2.sqrt() / (a_fro * x_norm + b_norm + 1e-30);
        worst_res = worst_res.max(res);
    }
    let pivot_ratio = if min_piv > 0.0 && min_piv.is_finite() {
        max_piv / min_piv
    } else {
        f64::INFINITY
    };
    Some((x_cols, worst_res, pivot_ratio))
}

fn solve_dense_checked_with_reg(
    label: &str,
    mat: &[C],
    rhs_cols: &[Vec<C>],
    p: usize,
    rep: &mut NetworkConversionReport,
) -> Vec<Vec<C>> {
    if let Some((sol, res, ratio)) = solve_dense_complex_multi_checked(mat, rhs_cols, p) {
        rep.max_solve_residual = rep.max_solve_residual.max(res);
        rep.max_pivot_ratio = rep.max_pivot_ratio.max(ratio);
        if ratio > 1e8 || res > 1e-8 {
            rep.warning_messages.push(format!(
                "ill-conditioned {label} network conversion (pivot ratio = {ratio:.2e}, residual = {res:.2e})"
            ));
        }
        return sol;
    }
    rep.used_regularization = true;
    let mut reg = mat.to_vec();
    for d in 0..p {
        reg[d * p + d] += c(1e-14, 0.0);
    }
    if let Some((sol, res, ratio)) = solve_dense_complex_multi_checked(&reg, rhs_cols, p) {
        rep.max_solve_residual = rep.max_solve_residual.max(res);
        rep.max_pivot_ratio = rep.max_pivot_ratio.max(ratio);
        rep.warning_messages.push(format!(
            "{label} network conversion encountered a singular matrix and used 1e-14 diagonal regularization (pivot ratio = {ratio:.2e})"
        ));
        sol
    } else {
        rep.singular_fallback = true;
        rep.max_solve_residual = 1.0;
        rep.max_pivot_ratio = f64::INFINITY;
        rep.warning_messages.push(format!(
            "{label} network conversion failed on a singular matrix and fell back to zeros"
        ));
        vec![vec![c(0.0, 0.0); p]; p]
    }
}

/// Computes `(passivity_max_sigma, reciprocity_rel_err)` of a $P \times P$ scattering matrix $\mathbf{S}$:
/// $\sigma_{\max}(\mathbf{S})$ and $\lVert \mathbf{S} - \mathbf{S}^T \rVert_F / \max(\lVert \mathbf{S} \rVert_F, 1.0)$.
pub fn s_matrix_passivity_and_reciprocity(s_mat: &[Vec<C>]) -> (f64, f64) {
    let p = s_mat.len();
    if p == 0 {
        return (0.0, 0.0);
    }
    let mut s_fro2 = 0.0_f64;
    let mut asym2 = 0.0_f64;
    let mut col_data = Vec::with_capacity(p * p);
    for col in 0..p {
        for r in 0..p {
            let s_rc = s_mat[r][col];
            col_data.push(s_rc);
            s_fro2 += s_rc.norm_sqr();
            asym2 += (s_rc - s_mat[col][r]).norm_sqr();
        }
    }
    let reciprocity = asym2.sqrt() / s_fro2.sqrt().max(1.0);
    let sigma_max = if p == 1 {
        s_mat[0][0].norm()
    } else if !col_data.iter().all(|v| v.is_finite()) {
        f64::INFINITY
    } else {
        crate::linalg::faer_glue::thin_svd_c64(&col_data, p, p)
            .ok()
            .and_then(|(s, _)| s.first().copied())
            .unwrap_or_else(|| s_fro2.sqrt())
    };
    (sigma_max, reciprocity)
}

/// Checked variant of [`s_to_z_and_y`] returning conditioning and singularity diagnostics.
pub fn s_to_z_and_y_checked(
    s_mat: &[Vec<C>],
    z0_ports: &[f64],
) -> (PortMatrixPairs, PortMatrixPairs, PortMatrixPairs, NetworkConversionReport) {
    let p = s_mat.len();
    let sqrt_z0: Vec<f64> = z0_ports.iter().map(|&z| z.max(1e-30).sqrt()).collect();
    let inv_sqrt_z0: Vec<f64> = sqrt_z0.iter().map(|&sz| 1.0 / sz).collect();

    let mut v_t = vec![c(0.0, 0.0); p * p];
    let mut i_t = vec![c(0.0, 0.0); p * p];
    let mut v_cols_t = vec![vec![c(0.0, 0.0); p]; p];
    let mut i_cols_t = vec![vec![c(0.0, 0.0); p]; p];

    for r in 0..p {
        for col in 0..p {
            let delta = if r == col { c(1.0, 0.0) } else { c(0.0, 0.0) };
            let v_rc = (delta + s_mat[r][col]) * sqrt_z0[r];
            let i_rc = (delta - s_mat[r][col]) * inv_sqrt_z0[r];
            v_t[col * p + r] = v_rc;
            i_t[col * p + r] = i_rc;
            v_cols_t[r][col] = v_rc;
            i_cols_t[r][col] = i_rc;
        }
    }

    let mut rep = NetworkConversionReport::default();
    let z_rows = solve_dense_checked_with_reg("S->Z", &i_t, &v_cols_t, p, &mut rep);
    let y_rows = solve_dense_checked_with_reg("S->Y", &v_t, &i_cols_t, p, &mut rep);

    let mut s_out = vec![vec![[0.0; 2]; p]; p];
    let mut z_out = vec![vec![[0.0; 2]; p]; p];
    let mut y_out = vec![vec![[0.0; 2]; p]; p];
    for r in 0..p {
        for col in 0..p {
            s_out[r][col] = [s_mat[r][col].re, s_mat[r][col].im];
            z_out[r][col] = [z_rows[r][col].re, z_rows[r][col].im];
            y_out[r][col] = [y_rows[r][col].re, y_rows[r][col].im];
        }
    }
    (s_out, z_out, y_out, rep)
}

/// Converts a $P \times P$ admittance matrix $\mathbf{Y}$ and positive reference impedances $Z_{0,p}$
/// into `(S, Z, Y)` arrays stored as `[re, im]` pairs, using
/// $\mathbf{S} = (\mathbf{I} + \mathbf{Y}_{\text{norm}})^{-1}(\mathbf{I} - \mathbf{Y}_{\text{norm}})$
/// with $(\mathbf{Y}_{\text{norm}})_{pq} = \sqrt{Z_{0,p}} Y_{pq} \sqrt{Z_{0,q}}$, and $\mathbf{Z} = \mathbf{Y}^{-1}$.
pub fn y_to_z_and_s(
    y_mat: &[Vec<C>],
    z0_ports: &[f64],
) -> (PortMatrixPairs, PortMatrixPairs, PortMatrixPairs) {
    let (s, z, y, _rep) = y_to_z_and_s_checked(y_mat, z0_ports);
    (s, z, y)
}

/// Checked variant of [`y_to_z_and_s`] returning conditioning and singularity diagnostics.
pub fn y_to_z_and_s_checked(
    y_mat: &[Vec<C>],
    z0_ports: &[f64],
) -> (PortMatrixPairs, PortMatrixPairs, PortMatrixPairs, NetworkConversionReport) {
    let p = y_mat.len();
    let sqrt_z0: Vec<f64> = z0_ports.iter().map(|&z| z.max(1e-30).sqrt()).collect();

    let mut y_dense = vec![c(0.0, 0.0); p * p];
    let mut ip_ynorm = vec![c(0.0, 0.0); p * p];
    let mut im_ynorm_cols = vec![vec![c(0.0, 0.0); p]; p];
    let mut eye_cols = vec![vec![c(0.0, 0.0); p]; p];

    for r in 0..p {
        eye_cols[r][r] = c(1.0, 0.0);
        for col in 0..p {
            y_dense[r * p + col] = y_mat[r][col];
            let yn = y_mat[r][col] * (sqrt_z0[r] * sqrt_z0[col]);
            let delta = if r == col { c(1.0, 0.0) } else { c(0.0, 0.0) };
            ip_ynorm[r * p + col] = delta + yn;
            im_ynorm_cols[col][r] = delta - yn;
        }
    }

    let mut rep = NetworkConversionReport::default();
    let s_cols = solve_dense_checked_with_reg("Y->S", &ip_ynorm, &im_ynorm_cols, p, &mut rep);
    let z_cols = solve_dense_checked_with_reg("Y->Z", &y_dense, &eye_cols, p, &mut rep);

    let mut s_out = vec![vec![[0.0; 2]; p]; p];
    let mut z_out = vec![vec![[0.0; 2]; p]; p];
    let mut y_out = vec![vec![[0.0; 2]; p]; p];
    for r in 0..p {
        for col in 0..p {
            let sij = s_cols[col][r];
            let zij = z_cols[col][r];
            let yij = y_mat[r][col];
            s_out[r][col] = [sij.re, sij.im];
            z_out[r][col] = [zij.re, zij.im];
            y_out[r][col] = [yij.re, yij.im];
        }
    }
    (s_out, z_out, y_out, rep)
}

/// Partitions `freqs_hz` into `[start..end)` sub-bands with `f_end / f_start <= 3.0` (`~2.5x` per band when `f_max / f_min > 3.0`).
#[allow(clippy::single_range_in_vec_init)]
pub fn partition_prom_subbands(freqs_hz: &[f64]) -> Vec<Range<usize>> {
    let n = freqs_hz.len();
    if n == 0 {
        return Vec::new();
    }
    if freqs_hz[n - 1] <= 3.0 * freqs_hz[0] {
        return vec![0..n];
    }
    let mut bands = Vec::new();
    let mut start = 0usize;
    while start < n {
        let f_cap = 2.5 * freqs_hz[start];
        let mut end = start + 1;
        while end < n && freqs_hz[end] <= f_cap {
            end += 1;
        }
        // Absorb trailing sliver if remaining span is <= 3.0x so the final sub-band is never degenerate.
        if end < n && freqs_hz[n - 1] <= 3.0 * freqs_hz[start] {
            end = n;
        }
        bands.push(start..end);
        start = end;
    }
    bands
}

/// Linear solver backend reusable across Hermite anchor derivatives (`Y`, `Y'`, `Y''`) and symbolic patterns.
pub trait PromLinearSolver {
    fn update_and_solve(
        &mut self,
        omega: f64,
        a_mat: &CCsr,
        rhs_cols: &[Vec<C>],
    ) -> Result<Vec<Vec<C>>, String>;

    fn solve_same_factor(&self, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String>;
}

/// Shared sparse direct / LBLT / BLR LU solver wrapper for both `nyx-ems` and `nyx-sms` driven solves.
pub enum SparseDrivenSolver {
    Direct {
        solver: Option<Box<DirectSolver<C>>>,
        requested: String,
        force_general_lu: bool,
    },
    Blr {
        solver: Option<Box<ShiftedBlrLu<C>>>,
        a_cur: CCsr,
        tol: f64,
    },
}

impl SparseDrivenSolver {
    pub fn new(requested: &str, blr_tol: f64) -> Self {
        Self::new_with_symmetry(requested, blr_tol, false)
    }

    /// Creates a driven sparse solver, routing to general `LU` when `is_non_symmetric` is `true` (`FR-34`).
    pub fn new_with_symmetry(requested: &str, blr_tol: f64, is_non_symmetric: bool) -> Self {
        if requested == "blr" {
            Self::Blr { solver: None, a_cur: CCsr::zeros(0, 0), tol: blr_tol }
        } else {
            Self::Direct {
                solver: None,
                requested: requested.to_string(),
                force_general_lu: is_non_symmetric,
            }
        }
    }

    pub fn update_matrix(&mut self, a_mat: &CCsr) -> Result<(), String> {
        match self {
            Self::Direct { solver, force_general_lu, .. } => match solver.as_mut() {
                Some(s) => s.update_matrix(a_mat),
                None => {
                    let ds = if *force_general_lu {
                        DirectSolver::new_general_lu(a_mat)?
                    } else {
                        DirectSolver::new(a_mat, false)?
                    };
                    *solver = Some(Box::new(ds));
                    Ok(())
                }
            },
            Self::Blr { solver, a_cur, tol } => match solver.as_mut() {
                Some(s) => {
                    *a_cur = CCsr::zeros(0, 0);
                    s.replace_matrices(a_mat, None).map_err(|e| format!("BLR matrix update failed: {e}"))?;
                    crate::par::on_factor_pool(|| s.set_shift(c(0.0, 0.0)))
                        .map_err(|e| format!("BLR numeric factorisation failed: {e}"))?;
                    *a_cur = a_mat.clone();
                    Ok(())
                }
                None => {
                    let opts = BlrOptions { tol: *tol, ..BlrOptions::default() };
                    let mut blr = ShiftedBlrLu::new(a_mat, None, &opts)
                        .map_err(|e| format!("BLR symbolic analysis failed: {e}"))?;
                    crate::par::on_factor_pool(|| blr.set_shift(c(0.0, 0.0)))
                        .map_err(|e| format!("BLR numeric factorisation failed: {e}"))?;
                    *solver = Some(Box::new(blr));
                    *a_cur = a_mat.clone();
                    Ok(())
                }
            },
        }
    }

    pub fn solve_many(&self, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String> {
        match self {
            Self::Direct { solver, .. } => {
                let s = solver.as_ref().ok_or("Direct solver has not been factorized yet")?;
                Ok(s.solve_many(rhs_cols))
            }
            Self::Blr { solver, a_cur, tol } => {
                let s = solver.as_ref().ok_or("BLR solver has not been factorized yet")?;
                if *tol <= 0.0 {
                    return Ok(s.solve_many(rhs_cols));
                }
                let n = s.dim();
                let mut out = Vec::with_capacity(rhs_cols.len());
                for b in rhs_cols {
                    let mut x = vec![c(0.0, 0.0); n];
                    s.solve_refined(a_cur, b, &mut x, 4, 1e-11);
                    out.push(x);
                }
                Ok(out)
            }
        }
    }

    /// Multi-RHS solve checking the relative backward residual against `a_mat` and dynamic `rcond`.
    pub fn solve_many_checked(&self, a_mat: &CCsr, rhs_cols: &[Vec<C>]) -> Result<(Vec<Vec<C>>, DirectSolveCheck), String> {
        match self {
            Self::Direct { solver, .. } => {
                let s = solver.as_ref().ok_or("Direct solver has not been factorized yet")?;
                Ok(s.solve_many_checked(a_mat, rhs_cols))
            }
            Self::Blr { solver, .. } => {
                let s = solver.as_ref().ok_or("BLR solver has not been factorized yet")?;
                let xs = self.solve_many(rhs_cols)?;
                let algebraic_residual = backward_residual_many(a_mat, rhs_cols, &xs);
                Ok((
                    xs,
                    DirectSolveCheck {
                        algebraic_residual,
                        rcond_estimate: None,
                        stats: s.stats(),
                        diagnostics: Default::default(),
                        cholesky_fallbacks: 0,
                        pivot_growth: 1.0,
                        solver_kind: "blr_lu",
                    },
                ))
            }
        }
    }

    pub fn symbolic_count(&self) -> usize {
        match self {
            Self::Direct { solver, .. } => solver.as_ref().map_or(0, |s| s.stats().symbolic_count),
            Self::Blr { solver, .. } => solver.as_ref().map_or(0, |s| s.stats().symbolic_count),
        }
    }

    pub fn backend_name(&self) -> String {
        match self {
            Self::Direct { solver, requested, .. } => {
                let Some(s) = solver.as_ref() else {
                    return "None".to_string();
                };
                let b = format_factor_backend(s.stats().backend);
                if requested == "lblt"
                    && !matches!(s.stats().backend, FactorBackend::LbltNd | FactorBackend::LbltAmd)
                {
                    format!("{b} (fallback from lblt)")
                } else {
                    b.to_string()
                }
            }
            Self::Blr { solver, .. } => solver
                .as_ref()
                .map_or_else(|| "None".to_string(), |s| format_factor_backend(s.stats().backend).to_string()),
        }
    }
}

impl PromLinearSolver for SparseDrivenSolver {
    fn update_and_solve(
        &mut self,
        _omega: f64,
        a_mat: &CCsr,
        rhs_cols: &[Vec<C>],
    ) -> Result<Vec<Vec<C>>, String> {
        self.update_matrix(a_mat)?;
        self.solve_many(rhs_cols)
    }

    fn solve_same_factor(&self, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String> {
        self.solve_many(rhs_cols)
    }
}

/// Affine frequency-domain operator $\hat{\mathbf{A}}(\omega) = \sum_m \alpha_m(\omega) \hat{\mathbf{A}}_m$
/// consumed by [`ImageProm`].
pub trait PromOperator: Sync {
    fn n_dofs(&self) -> usize;
    fn n_ports(&self) -> usize;
    fn terms_hat(&self) -> &[CCsr];
    fn term_coeffs(&self, omega: f64) -> Vec<(C, C, C)>;
    fn d_ruiz(&self) -> &[f64];
    fn port_b_hat(&self, port: usize) -> Vec<C>;
    /// Optional non-symmetric port functional vector for response recovery (`b_recv != b_drive`).
    fn port_b_recv_hat(&self, _port: usize) -> Option<Vec<C>> {
        None
    }
    /// Scalar frequency multiplier $(\gamma_q(\omega), \gamma_q'(\omega), \gamma_q''(\omega))$
    /// such that $\text{RHS}_q(\omega) = \gamma_q(\omega) \hat{\mathbf{b}}_q$.
    fn port_gamma(&self, port: usize, omega: f64) -> (C, C, C);
    /// Whether $\hat{\mathbf{A}}_m = \hat{\mathbf{A}}_m^T$ (complex-symmetric) so only upper-triangle inner products are needed.
    fn is_symmetric(&self) -> bool {
        true
    }
    /// Optional disjoint index ranges partitioning snapshots before SVD truncation (e.g. `(u, phi)` blocks).
    fn block_splits(&self) -> Option<&[Range<usize>]> {
        None
    }

    fn assemble_matrix_at(&self, omega: f64) -> CCsr {
        let coeffs = self.term_coeffs(omega);
        let terms = self.terms_hat();
        let mut mat = terms[0].scaled(coeffs[0].0);
        for (m, &(alpha, _, _)) in coeffs.iter().enumerate().skip(1) {
            if terms[m].nnz() > 0 {
                mat = mat.add_scaled(&terms[m], alpha);
            }
        }
        mat
    }

    fn apply_matrix_deriv(&self, omega: f64, order_deriv: usize, x: &[C]) -> Vec<C> {
        let coeffs = self.term_coeffs(omega);
        let terms = self.terms_hat();
        let n = self.n_dofs();
        let mut out = vec![c(0.0, 0.0); n];
        for (m, &(_, d1, d2)) in coeffs.iter().enumerate() {
            let coeff = if order_deriv == 1 { d1 } else { d2 };
            if (coeff.re == 0.0 && coeff.im == 0.0) || terms[m].nnz() == 0 {
                continue;
            }
            let mx = terms[m].matvec(x);
            for (yi, &mxi) in out.iter_mut().zip(&mx) {
                *yi += coeff * mxi;
            }
        }
        out
    }
}

fn transpose_ccsr(a: &CCsr) -> CCsr {
    let mut coo = CCoo::new(a.ncols(), a.nrows());
    for i in 0..a.nrows() {
        let (cols, vals) = a.row(i);
        for (&j, &v) in cols.iter().zip(vals) {
            coo.push(j as usize, i, v);
        }
    }
    coo.to_csr()
}

/// Real Hermite Tangential Galerkin PROM with isometric QR image-space basis $\mathbf{U}$ for exact
/// microsecond full-space residual certification across all sweep frequencies.
pub struct ImageProm<'a, O: PromOperator> {
    pub op: &'a O,
    pub q_cols: Vec<Vec<f64>>,
    pub a_r: Vec<Vec<C>>,
    pub b_r_drive: Vec<Vec<C>>,
    pub b_r_recv: Vec<Vec<C>>,
    pub has_distinct_recv: bool,
    pub u_cols: Vec<Vec<f64>>,
    pub u_b: Vec<Vec<C>>,
    pub b_norm2: Vec<f64>,
    pub u_aq: Vec<Vec<Vec<C>>>,
    terms_hat_t: Vec<CCsr>,
    b_drive_hat: Vec<Vec<C>>,
    b_recv_hat: Vec<Vec<C>>,
}

impl<'a, O: PromOperator> ImageProm<'a, O> {
    pub fn new(op: &'a O) -> Self {
        let terms = op.terms_hat();
        let m_terms = terms.len();
        let p_ports = op.n_ports();
        let b_drive_hat: Vec<Vec<C>> = (0..p_ports).map(|p| op.port_b_hat(p)).collect();
        let mut has_distinct_recv = false;
        let b_recv_hat: Vec<Vec<C>> = (0..p_ports)
            .map(|p| {
                if let Some(br) = op.port_b_recv_hat(p) {
                    has_distinct_recv = true;
                    br
                } else {
                    b_drive_hat[p].clone()
                }
            })
            .collect();

        let mut u_cols: Vec<Vec<f64>> = Vec::new();
        for bh in &b_drive_hat {
            let re_v: Vec<f64> = bh.iter().map(|z| z.re).collect();
            let im_v: Vec<f64> = bh.iter().map(|z| z.im).collect();
            append_l2_orthonormal(&mut u_cols, &re_v);
            append_l2_orthonormal(&mut u_cols, &im_v);
        }
        let s = u_cols.len();
        let mut u_b = Vec::with_capacity(p_ports);
        let mut b_norm2 = Vec::with_capacity(p_ports);
        for bh in &b_drive_hat {
            b_norm2.push(bh.iter().map(|z| z.norm_sqr()).sum());
            let coords: Vec<C> = u_cols.iter().map(|u| dot_real_complex(u, bh)).collect();
            debug_assert_eq!(coords.len(), s);
            u_b.push(coords);
        }
        let terms_hat_t = if op.is_symmetric() { Vec::new() } else { terms.iter().map(transpose_ccsr).collect() };

        Self {
            op,
            q_cols: Vec::new(),
            a_r: vec![Vec::new(); m_terms],
            b_r_drive: Vec::new(),
            b_r_recv: Vec::new(),
            has_distinct_recv,
            u_cols,
            u_b,
            b_norm2,
            u_aq: vec![Vec::new(); m_terms],
            terms_hat_t,
            b_drive_hat,
            b_recv_hat,
        }
    }

    #[inline]
    pub fn dim(&self) -> usize {
        self.q_cols.len()
    }

    /// Returns projected port vectors `b_r_recv[q_basis][port]` used for port recovery ($\mathbf{b}_p^T \hat{\mathbf{x}}_q$).
    #[inline]
    pub fn b_r_recv(&self) -> &[Vec<C>] {
        if self.has_distinct_recv {
            &self.b_r_recv
        } else {
            &self.b_r_drive
        }
    }

    /// Appends a real candidate vector to $\mathbf{Q}$ (if not already in the span) and updates
    /// $\mathbf{A}_{m,r}$, $\mathbf{B}_r$, and the isometric image-space coordinates in $\mathbf{U}$.
    pub fn add_real_vector(&mut self, cand: &[f64]) -> bool {
        if !append_l2_orthonormal(&mut self.q_cols, cand) {
            return false;
        }
        let r_old = self.q_cols.len() - 1;
        let r_new = r_old + 1;
        let q_new = &self.q_cols[r_old];
        let terms = self.op.terms_hat();
        let m_terms = terms.len();
        let nred = self.op.n_dofs();
        let is_sym = self.op.is_symmetric();

        let b_row_drive: Vec<C> = self.b_drive_hat.iter().map(|bh| dot_real_complex(q_new, bh)).collect();
        self.b_r_drive.push(b_row_drive);
        if self.has_distinct_recv {
            let b_row_recv: Vec<C> = self.b_recv_hat.iter().map(|bh| dot_real_complex(q_new, bh)).collect();
            self.b_r_recv.push(b_row_recv);
        }

        let q_new_c: Vec<C> = q_new.iter().map(|&x| c(x, 0.0)).collect();
        for m in 0..m_terms {
            let ym = if terms[m].nnz() > 0 { terms[m].matvec(&q_new_c) } else { vec![c(0.0, 0.0); nred] };

            let mut grown = vec![c(0.0, 0.0); r_new * r_new];
            for i in 0..r_old {
                for j in 0..r_old {
                    grown[i * r_new + j] = self.a_r[m][i * r_old + j];
                }
            }
            for i in 0..r_new {
                let qi = &self.q_cols[i];
                let dot_col = dot_real_complex(qi, &ym);
                grown[i * r_new + r_old] = dot_col;
                if is_sym {
                    grown[r_old * r_new + i] = dot_col;
                }
            }
            if !is_sym {
                let ym_t = if self.terms_hat_t[m].nnz() > 0 {
                    self.terms_hat_t[m].matvec(&q_new_c)
                } else {
                    vec![c(0.0, 0.0); nred]
                };
                for j in 0..r_new {
                    grown[r_old * r_new + j] = dot_real_complex(&self.q_cols[j], &ym_t);
                }
            }
            self.a_r[m] = grown;

            let re_y: Vec<f64> = ym.iter().map(|z| z.re).collect();
            let im_y: Vec<f64> = ym.iter().map(|z| z.im).collect();
            append_l2_orthonormal(&mut self.u_cols, &re_y);
            append_l2_orthonormal(&mut self.u_cols, &im_y);

            let s_now = self.u_cols.len();
            let mut coords = Vec::with_capacity(s_now);
            for u in &self.u_cols {
                coords.push(dot_real_complex(u, &ym));
            }
            self.u_aq[m].push(coords);
        }

        let s_now = self.u_cols.len();
        for ub in &mut self.u_b {
            ub.resize(s_now, c(0.0, 0.0));
        }
        for m in 0..m_terms {
            for col in &mut self.u_aq[m] {
                col.resize(s_now, c(0.0, 0.0));
            }
        }
        true
    }

    /// Splits a complex snapshot vector `z` into phase-aligned `Re(z)` and `Im(z)` and appends both to $\mathbf{Q}$.
    pub fn add_complex_snapshot(&mut self, z: &[C]) {
        let big = z.iter().copied().fold(c(0.0, 0.0), |m, v| if v.norm_sqr() > m.norm_sqr() { v } else { m });
        if big.norm_sqr() <= 1e-300 {
            return;
        }
        let rot = big.conj() / big.norm();
        let re_v: Vec<f64> = z.iter().map(|&v| (v * rot).re).collect();
        let im_v: Vec<f64> = z.iter().map(|&v| (v * rot).im).collect();
        self.add_real_vector(&re_v);
        self.add_real_vector(&im_v);
    }

    fn ingest_anchor_block(&mut self, cols: &[Vec<C>], range: Option<Range<usize>>) -> Result<(), String> {
        let n = self.op.n_dofs();
        let slice_cols: Vec<Vec<C>> = cols
            .iter()
            .filter_map(|col| {
                let seg = if let Some(ref rg) = range { col[rg.clone()].to_vec() } else { col.clone() };
                let nrm = seg.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
                (nrm > 1e-30 && nrm.is_finite()).then(|| {
                    let inv = 1.0 / nrm;
                    seg.into_iter().map(|v| v * inv).collect()
                })
            })
            .collect();
        if slice_cols.is_empty() {
            return Ok(());
        }
        let (s_vals, u_cols) = thin_svd_c64(&slice_cols)?;
        if let Some(&s0) = s_vals.first() {
            if s0 > 0.0 {
                let cutoff = 1e-6 * s0;
                for (&sk, uk) in s_vals.iter().zip(&u_cols) {
                    if sk >= cutoff {
                        if let Some(ref rg) = range {
                            let mut padded = vec![c(0.0, 0.0); n];
                            padded[rg.clone()].copy_from_slice(uk);
                            self.add_complex_snapshot(&padded);
                        } else {
                            self.add_complex_snapshot(uk);
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Performs one full HDM numeric factorization at angular frequency `omega_a` (rad/s),
    /// solves for `Y`, `Y'`, and `Y''` on the same factorization, compresses via thin SVD, and expands $\mathbf{Q}$.
    pub fn add_hermite_anchor<S: PromLinearSolver>(&mut self, omega_a: f64, solver: &mut S) -> Result<(), String> {
        let a_mat = self.op.assemble_matrix_at(omega_a);
        let p_ports = self.op.n_ports();
        let nred = self.op.n_dofs();

        let mut rhs_0 = Vec::with_capacity(p_ports);
        let mut gammas = Vec::with_capacity(p_ports);
        for q in 0..p_ports {
            let g = self.op.port_gamma(q, omega_a);
            gammas.push(g);
            rhs_0.push(self.b_drive_hat[q].iter().map(|&b| g.0 * b).collect::<Vec<C>>());
        }

        // 1. Value solve: A(omega_a) Y = RHS(omega_a)
        let y0_cols = solver.update_and_solve(omega_a, &a_mat, &rhs_0)?;

        // 2. First derivative solve on the SAME factorization:
        //    A(omega_a) Y' = RHS'(omega_a) - A'(omega_a) Y
        let mut rhs_1 = Vec::with_capacity(p_ports);
        for q in 0..p_ports {
            let a1_y0 = self.op.apply_matrix_deriv(omega_a, 1, &y0_cols[q]);
            let g1 = gammas[q].1;
            let b_hat = &self.b_drive_hat[q];
            let col: Vec<C> = (0..nred).map(|i| g1 * b_hat[i] - a1_y0[i]).collect();
            rhs_1.push(col);
        }
        let y1_cols = solver.solve_same_factor(&rhs_1)?;

        // 3. Second derivative solve on the SAME factorization:
        //    A(omega_a) Y'' = RHS''(omega_a) - 2 A'(omega_a) Y' - A''(omega_a) Y
        let mut rhs_2 = Vec::with_capacity(p_ports);
        for q in 0..p_ports {
            let a1_y1 = self.op.apply_matrix_deriv(omega_a, 1, &y1_cols[q]);
            let a2_y0 = self.op.apply_matrix_deriv(omega_a, 2, &y0_cols[q]);
            let g2 = gammas[q].2;
            let b_hat = &self.b_drive_hat[q];
            let col: Vec<C> = (0..nred).map(|i| g2 * b_hat[i] - c(2.0, 0.0) * a1_y1[i] - a2_y0[i]).collect();
            rhs_2.push(col);
        }
        let y2_cols = solver.solve_same_factor(&rhs_2)?;

        let all_cols: Vec<Vec<C>> = y0_cols.into_iter().chain(y1_cols).chain(y2_cols).collect();
        if let Some(splits) = self.op.block_splits() {
            for rg in splits {
                self.ingest_anchor_block(&all_cols, Some(rg.clone()))?;
            }
        } else {
            self.ingest_anchor_block(&all_cols, None)?;
        }
        Ok(())
    }

    /// Solves the $r \times r$ reduced system at `omega` (rad/s) for all $P$ ports and returns
    /// `(x_r_cols, exact_relative_residual)`.
    pub fn solve_reduced_and_residual(&self, omega: f64) -> Option<(Vec<Vec<C>>, f64)> {
        let r = self.dim();
        let p_ports = self.op.n_ports();
        let m_terms = self.op.terms_hat().len();
        let coeffs = self.op.term_coeffs(omega);

        let mut t_r = vec![c(0.0, 0.0); r * r];
        for m in 0..m_terms {
            let alpha = coeffs[m].0;
            if alpha.re == 0.0 && alpha.im == 0.0 {
                continue;
            }
            let am = &self.a_r[m];
            for idx in 0..(r * r) {
                t_r[idx] += alpha * am[idx];
            }
        }

        let mut gammas = Vec::with_capacity(p_ports);
        let mut rhs_r = Vec::with_capacity(p_ports);
        for q in 0..p_ports {
            let (g0, _, _) = self.op.port_gamma(q, omega);
            gammas.push(g0);
            let col: Vec<C> = (0..r).map(|i| g0 * self.b_r_drive[i][q]).collect();
            rhs_r.push(col);
        }

        let xr_cols = solve_dense_complex_multi(&t_r, &rhs_r, r)?;

        let s = self.u_cols.len();
        let mut res_num2 = 0.0;
        let mut res_den2 = 0.0;
        let mut r_vec = vec![c(0.0, 0.0); s];

        for q in 0..p_ports {
            let g0 = gammas[q];
            res_den2 += g0.norm_sqr() * self.b_norm2[q];
            for (k, &ub_k) in self.u_b[q].iter().enumerate() {
                r_vec[k] = -g0 * ub_k;
            }
            let xrq = &xr_cols[q];
            for m in 0..m_terms {
                let alpha = coeffs[m].0;
                if alpha.re == 0.0 && alpha.im == 0.0 {
                    continue;
                }
                let u_aq_m = &self.u_aq[m];
                for i in 0..r {
                    let weight = alpha * xrq[i];
                    if weight.re == 0.0 && weight.im == 0.0 {
                        continue;
                    }
                    let col_u = &u_aq_m[i];
                    for k in 0..s {
                        r_vec[k] += weight * col_u[k];
                    }
                }
            }
            for &rv in &r_vec {
                res_num2 += rv.norm_sqr();
            }
        }

        let rel_res = (res_num2 / res_den2.max(1e-300)).sqrt();
        Some((xr_cols, rel_res))
    }

    /// Reconstructs the unscaled reduced-DOF physical solution vector $\mathbf{x} = \mathbf{D} \mathbf{Q} \mathbf{x}_r$.
    pub fn reconstruct_unscaled(&self, xr: &[C]) -> Vec<C> {
        let nred = self.op.n_dofs();
        let mut y = vec![c(0.0, 0.0); nred];
        for (i, &c_val) in xr.iter().enumerate() {
            let qi = &self.q_cols[i];
            for k in 0..nred {
                y[k] += c_val * qi[k];
            }
        }
        for (yk, &dk) in y.iter_mut().zip(self.op.d_ruiz()) {
            *yk *= dk;
        }
        y
    }
}

/// Summary metrics produced by [`run_greedy_prom_sweep`].
pub struct PromSweepSummary {
    pub hdm_solves: usize,
    pub rom_dim: usize,
    pub anchor_freqs_hz: Vec<f64>,
    pub max_residual: f64,
    pub prom_residuals: Vec<f64>,
    pub rom_sweep_seconds: f64,
}

/// Runs a multi-band adaptive Hermite/Galerkin PROM sweep across `freqs_hz`, certifying the full-space
/// residual at every point in `U` and calling `on_band` with each converged band's `ImageProm` and reduced
/// solutions `xr_band[k_local][port][mode]`.
#[allow(clippy::too_many_arguments)]
pub fn run_greedy_prom_sweep<O, S, F>(
    op: &O,
    solver: &mut S,
    freqs_hz: &[f64],
    prom_tol: f64,
    max_hdm_solves: usize,
    strict_prom_tol: bool,
    seed_freqs_hz: &[f64],
    mut on_band: F,
) -> Result<PromSweepSummary, String>
where
    O: PromOperator,
    S: PromLinearSolver,
    F: FnMut(Range<usize>, &[f64], &ImageProm<'_, O>, &[Vec<Vec<C>>]) -> Result<(), String>,
{
    let omegas: Vec<f64> = freqs_hz.iter().map(|&f| 2.0 * std::f64::consts::PI * f).collect();
    let bands = partition_prom_subbands(freqs_hz);
    let n_bands = bands.len().max(1);
    let total_max_solves = max_hdm_solves.max(1);
    let p_ports = op.n_ports();
    const SWEEP_CHUNK: usize = 32;

    let mut hdm_solves = 0usize;
    let mut rom_dim = 0usize;
    let mut anchor_freqs_hz = Vec::new();
    let mut max_residual = 0.0_f64;
    let mut prom_residuals = vec![0.0_f64; freqs_hz.len()];
    let mut rom_sweep_seconds = 0.0_f64;

    for (bi, band) in bands.into_iter().enumerate() {
        let b_freqs = &freqs_hz[band.clone()];
        let b_omegas = &omegas[band.clone()];
        let nb_pts = b_freqs.len();
        let remaining = total_max_solves.saturating_sub(hdm_solves).max(1);
        let bands_left = (n_bands - bi).max(1);
        let max_solves_band =
            if n_bands == 1 { total_max_solves } else { remaining.div_ceil(bands_left).max(3.min(remaining)) };

        let mut prom = ImageProm::new(op);
        let mut band_solves = 0usize;

        let mut base_indices: Vec<usize> =
            if nb_pts == 1 || (b_freqs[nb_pts - 1] - b_freqs[0]).abs() <= 1e-12 * b_freqs[0] {
                vec![0]
            } else if nb_pts == 2 {
                vec![0, 1]
            } else {
                vec![0, (nb_pts - 1) / 2, nb_pts - 1]
            };
        base_indices.dedup();

        let mut initial_anchors: Vec<(Option<usize>, f64)> =
            base_indices.iter().map(|&idx| (Some(idx), b_freqs[idx])).collect();
        for &sf in seed_freqs_hz {
            if n_bands > 1 && (sf < b_freqs[0] - 1e-9 || sf > b_freqs[nb_pts - 1] + 1e-9) {
                continue;
            }
            let mut best_k = 0usize;
            let mut best_diff = f64::INFINITY;
            for (k, &fk) in b_freqs.iter().enumerate() {
                let diff = (fk - sf).abs();
                if diff < best_diff {
                    best_diff = diff;
                    best_k = k;
                }
            }
            if best_diff <= 1e-9 * sf.max(1.0) {
                if !initial_anchors.iter().any(|&(oi, _)| oi == Some(best_k)) {
                    initial_anchors.push((Some(best_k), b_freqs[best_k]));
                }
            } else if !initial_anchors.iter().any(|&(_, f)| (f - sf).abs() <= 1e-9 * sf.max(1.0)) {
                initial_anchors.push((None, sf));
            }
        }

        let mut anchor_indices: Vec<usize> = Vec::new();
        for (opt_idx, f_anchor) in initial_anchors.into_iter().take(max_solves_band) {
            crate::progress::check()?;
            prom.add_hermite_anchor(2.0 * std::f64::consts::PI * f_anchor, solver)?;
            band_solves += 1;
            hdm_solves += 1;
            anchor_freqs_hz.push(f_anchor);
            if let Some(idx) = opt_idx {
                anchor_indices.push(idx);
            }
        }

        let mut xr_band: Vec<Vec<Vec<C>>>;
        let band_worst_res;
        let band_res: Vec<f64>;
        loop {
            crate::progress::check()?;
            let t_sweep = Instant::now();
            let dim_now = prom.dim();
            let prom_ref = &prom;
            let sweep_chunks = crate::par::map_chunks(b_omegas, SWEEP_CHUNK, |chunk| {
                chunk
                    .iter()
                    .map(|&w| match prom_ref.solve_reduced_and_residual(w) {
                        Some(pair) => pair,
                        None => (vec![vec![c(0.0, 0.0); dim_now]; p_ports], f64::INFINITY),
                    })
                    .collect::<Vec<_>>()
            });
            xr_band = Vec::with_capacity(nb_pts);
            let mut cur_res = Vec::with_capacity(nb_pts);
            let mut worst_res = 0.0_f64;
            let mut worst_cand_res = -1.0_f64;
            let mut worst_cand_idx = None;
            for (k, (xr_k, res_k)) in sweep_chunks.into_iter().flatten().enumerate() {
                if res_k > worst_res {
                    worst_res = res_k;
                }
                if !anchor_indices.contains(&k) && res_k > worst_cand_res {
                    worst_cand_res = res_k;
                    worst_cand_idx = Some(k);
                }
                xr_band.push(xr_k);
                cur_res.push(res_k);
            }
            rom_sweep_seconds += t_sweep.elapsed().as_secs_f64();

            if worst_res <= prom_tol || band_solves >= max_solves_band || hdm_solves >= total_max_solves {
                band_worst_res = worst_res;
                band_res = cur_res;
                break;
            }
            let Some(next_idx) = worst_cand_idx else {
                band_worst_res = worst_res;
                band_res = cur_res;
                break;
            };
            anchor_indices.push(next_idx);
            prom.add_hermite_anchor(b_omegas[next_idx], solver)?;
            band_solves += 1;
            hdm_solves += 1;
            anchor_freqs_hz.push(b_freqs[next_idx]);
        }

        max_residual = max_residual.max(band_worst_res);
        rom_dim += prom.dim();
        prom_residuals[band.clone()].copy_from_slice(&band_res);
        on_band(band, b_omegas, &prom, &xr_band)?;
    }

    if strict_prom_tol && max_residual > prom_tol {
        return Err(format!(
            "PROM did not converge to prom_tol = {prom_tol:.3e} within max_hdm_solves = {total_max_solves} (max_prom_residual = {max_residual:.3e})"
        ));
    }

    Ok(PromSweepSummary { hdm_solves, rom_dim, anchor_freqs_hz, max_residual, prom_residuals, rom_sweep_seconds })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::equilibrate::{ruiz_symmetric_terms, scale_symmetric, scale_vec, RuizTerm, DEFAULT_RUIZ_PASSES};

    struct SyntheticQuadOp {
        n: usize,
        p: usize,
        d_ruiz: Vec<f64>,
        terms_hat: Vec<CCsr>,
        b_hat: Vec<Vec<C>>,
    }

    impl PromOperator for SyntheticQuadOp {
        fn n_dofs(&self) -> usize {
            self.n
        }

        fn n_ports(&self) -> usize {
            self.p
        }

        fn terms_hat(&self) -> &[CCsr] {
            &self.terms_hat
        }

        fn term_coeffs(&self, omega: f64) -> Vec<(C, C, C)> {
            vec![
                (c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
                (c(0.0, omega), c(0.0, 1.0), c(0.0, 0.0)),
                (c(-omega * omega, 0.0), c(-2.0 * omega, 0.0), c(-2.0, 0.0)),
            ]
        }

        fn d_ruiz(&self) -> &[f64] {
            &self.d_ruiz
        }

        fn port_b_hat(&self, port: usize) -> Vec<C> {
            self.b_hat[port].clone()
        }

        fn port_gamma(&self, _port: usize, omega: f64) -> (C, C, C) {
            (c(1.0, 0.05 * omega / 1e9), c(0.0, 0.05 / 1e9), c(0.0, 0.0))
        }
    }

    #[test]
    fn hermite_prom_matches_hdm_to_1e_minus_14_and_certifies_residual() {
        // 1D-like tridiagonal mechanical+lossy chain of 60 DOFs with 2 ports and wide eigenvalue scale contrast.
        let n = 60usize;
        let p = 2usize;
        let mut coo_k = CCoo::new(n, n);
        let mut coo_c = CCoo::new(n, n);
        let mut coo_m = CCoo::new(n, n);
        for i in 0..n {
            let k_diag = 4.0e11 * (1.0 + 0.02 * (i as f64 / n as f64));
            let c_diag = 5.0 * (1.0 + 0.01 * i as f64);
            let m_diag = 2.0e-9 * (1.0 + 0.03 * ((n - i) as f64 / n as f64));
            coo_k.push(i, i, c(k_diag, 2e8));
            coo_c.push(i, i, c(c_diag, 0.0));
            coo_m.push(i, i, c(m_diag, 0.0));
            if i + 1 < n {
                coo_k.push(i, i + 1, c(-1.6e11, -5e7));
                coo_k.push(i + 1, i, c(-1.6e11, -5e7));
                coo_c.push(i, i + 1, c(-1.5, 0.0));
                coo_c.push(i + 1, i, c(-1.5, 0.0));
            }
        }
        let k = coo_k.to_csr();
        let cc = coo_c.to_csr();
        let m = coo_m.to_csr();
        let b0: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else if i == 1 { -0.3 } else { 0.0 }).collect();
        let b1_c: Vec<C> = (0..n)
            .map(|i| if i == n - 1 { c(1.0, 0.25) } else if i == n - 2 { c(-0.2, 0.1) } else { c(0.0, 0.0) })
            .collect();

        let f_mid = 1.0e9;
        let w_mid = 2.0 * std::f64::consts::PI * f_mid;
        let d_ruiz = ruiz_symmetric_terms(
            n,
            &[RuizTerm::Complex(&k, 1.0), RuizTerm::Complex(&cc, w_mid), RuizTerm::Complex(&m, w_mid * w_mid)],
            DEFAULT_RUIZ_PASSES,
        );
        let terms_hat = vec![scale_symmetric(&k, &d_ruiz), scale_symmetric(&cc, &d_ruiz), scale_symmetric(&m, &d_ruiz)];
        let b0_hat: Vec<C> = scale_vec(&b0, &d_ruiz).into_iter().map(|v| c(v, 0.0)).collect();
        let b1_hat: Vec<C> = b1_c.iter().zip(&d_ruiz).map(|(&bi, &di)| bi * di).collect();

        let op = SyntheticQuadOp { n, p, d_ruiz: d_ruiz.clone(), terms_hat, b_hat: vec![b0_hat, b1_hat] };
        let freqs: Vec<f64> = (0..25).map(|i| 0.8e9 + (0.4e9) * (i as f64) / 24.0).collect();

        // Check 1-anchor (large residual) QR residual accuracy vs full-space exact residual without cancellation floor.
        {
            let mut op1_solver = SparseDrivenSolver::new("lblt", 0.0);
            let mut prom1 = ImageProm::new(&op);
            prom1.add_hermite_anchor(2.0 * std::f64::consts::PI * freqs[0], &mut op1_solver).unwrap();
            let w_test = 2.0 * std::f64::consts::PI * freqs[freqs.len() - 1];
            let (xr_test, qr_res) = prom1.solve_reduced_and_residual(w_test).unwrap();
            let a_hat = op.assemble_matrix_at(w_test);
            let mut num2 = 0.0;
            let mut den2 = 0.0;
            for q in 0..p {
                let g0 = op.port_gamma(q, w_test).0;
                let rhs_q: Vec<C> = op.port_b_hat(q).into_iter().map(|b| b * g0).collect();
                let x_hat_q: Vec<C> =
                    prom1.reconstruct_unscaled(&xr_test[q]).iter().zip(&op.d_ruiz).map(|(&xi, &di)| xi / di).collect();
                let ax = a_hat.matvec(&x_hat_q);
                for i in 0..n {
                    num2 += (ax[i] - rhs_q[i]).norm_sqr();
                    den2 += rhs_q[i].norm_sqr();
                }
            }
            let exact_res = (num2 / den2).sqrt();
            assert!(
                (qr_res - exact_res).abs() <= 1e-12 * exact_res + 1e-14,
                "1-anchor QR residual mismatch: qr={qr_res:.3e} vs exact={exact_res:.3e}"
            );
        }

        let mut solver = SparseDrivenSolver::new("lblt", 0.0);
        let mut prom_x_by_freq = vec![Vec::new(); freqs.len()];
        let summary = run_greedy_prom_sweep(
            &op,
            &mut solver,
            &freqs,
            1e-13,
            12,
            false,
            &[],
            |band, b_omegas, prom, xr_band| {
                for (k_loc, &w) in b_omegas.iter().enumerate() {
                    let (_, qr_res) = prom.solve_reduced_and_residual(w).unwrap();
                    let a_hat = op.assemble_matrix_at(w);
                    let mut num2 = 0.0;
                    let mut den2 = 0.0;
                    for q in 0..p {
                        let g0 = op.port_gamma(q, w).0;
                        let rhs_q: Vec<C> = op.port_b_hat(q).into_iter().map(|b| b * g0).collect();
                        let x_hat_q: Vec<C> = prom
                            .reconstruct_unscaled(&xr_band[k_loc][q])
                            .iter()
                            .zip(&op.d_ruiz)
                            .map(|(&xi, &di)| xi / di)
                            .collect();
                        let ax = a_hat.matvec(&x_hat_q);
                        for i in 0..n {
                            num2 += (ax[i] - rhs_q[i]).norm_sqr();
                            den2 += rhs_q[i].norm_sqr();
                        }
                    }
                    let exact_res = (num2 / den2).sqrt();
                    assert!(
                        (qr_res - exact_res).abs() <= 1e-12 * exact_res + 1e-13,
                        "QR residual mismatch: qr={qr_res:.3e} vs exact={exact_res:.3e}"
                    );
                }
                for (k_loc, xr_k) in xr_band.iter().enumerate() {
                    prom_x_by_freq[band.start + k_loc] =
                        xr_k.iter().map(|xrq| prom.reconstruct_unscaled(xrq)).collect();
                }
                Ok(())
            },
        )
        .unwrap();
        assert!(summary.max_residual < 1e-13, "max_residual = {:.3e}", summary.max_residual);
        assert!(summary.hdm_solves <= 10, "converged in {} anchors", summary.hdm_solves);

        // Compare against exact full HDM solve at every sweep point.
        for (k_idx, &f) in freqs.iter().enumerate() {
            let w = 2.0 * std::f64::consts::PI * f;
            let a_hat = op.assemble_matrix_at(w);
            let hdm_solver = DirectSolver::new(&a_hat, false).unwrap();
            let rhs: Vec<Vec<C>> =
                (0..p).map(|q| op.port_b_hat(q).into_iter().map(|b| b * op.port_gamma(q, w).0).collect()).collect();
            let x_hat_cols = hdm_solver.solve_many(&rhs);
            for q in 0..p {
                let x_hdm: Vec<C> = x_hat_cols[q].iter().zip(&d_ruiz).map(|(&xi, &di)| xi * di).collect();
                let x_prom = &prom_x_by_freq[k_idx][q];
                let mut err2 = 0.0;
                let mut ref2 = 0.0;
                for i in 0..n {
                    err2 += (x_prom[i] - x_hdm[i]).norm_sqr();
                    ref2 += x_hdm[i].norm_sqr();
                }
                let rel_err = (err2 / ref2.max(1e-300)).sqrt();
                assert!(rel_err < 1e-12, "freq {f}: rel_err = {rel_err:.3e}");
            }
        }
    }

    #[test]
    fn y_and_s_conversions_round_trip_and_preserve_reciprocity() {
        let y = vec![
            vec![c(0.012, 0.035), c(-0.002, -0.011)],
            vec![c(-0.002, -0.011), c(0.015, 0.029)],
        ];
        let z0 = [50.0, 75.0];
        let (s_arr, z_arr, y_arr) = y_to_z_and_s(&y, &z0);
        let s_c: Vec<Vec<C>> = s_arr.iter().map(|r| r.iter().map(|&[re, im]| c(re, im)).collect()).collect();
        let (s2, z2, y2) = s_to_z_and_y(&s_c, &z0);

        for i in 0..2 {
            for j in 0..2 {
                assert!((s_arr[i][j][0] - s2[i][j][0]).abs() < 1e-14);
                assert!((s_arr[i][j][1] - s2[i][j][1]).abs() < 1e-14);
                assert!((z_arr[i][j][0] - z2[i][j][0]).abs() < 1e-12);
                assert!((z_arr[i][j][1] - z2[i][j][1]).abs() < 1e-12);
                assert!((y_arr[i][j][0] - y2[i][j][0]).abs() < 1e-13);
                assert!((y_arr[i][j][1] - y2[i][j][1]).abs() < 1e-13);
            }
        }
        // Reciprocity: Z_{01} == Z_{10}, S_{01}/sqrt(...) -> actually S^T == S for Z0-normalized S!
        assert!((z_arr[0][1][0] - z_arr[1][0][0]).abs() < 1e-13);
        assert!((z_arr[0][1][1] - z_arr[1][0][1]).abs() < 1e-13);
        assert!((s_arr[0][1][0] - s_arr[1][0][0]).abs() < 1e-14);
        assert!((s_arr[0][1][1] - s_arr[1][0][1]).abs() < 1e-14);
    }

    #[test]
    fn checked_port_conversions_and_passivity_report_ill_conditioning() {
        // S = [[1, 0], [0, -1]] makes (I - S) and (I + S) rank-deficient -> triggers explicit warnings.
        let s_open = vec![vec![c(1.0, 0.0), c(0.0, 0.0)], vec![c(0.0, 0.0), c(-1.0, 0.0)]];
        let (_s, _z, _y, rep) = s_to_z_and_y_checked(&s_open, &[50.0, 50.0]);
        assert!(
            rep.used_regularization || rep.singular_fallback || !rep.warning_messages.is_empty(),
            "expected singular S->Z/Y warning report, got {rep:?}"
        );

        // Non-passive and non-reciprocal 2x2 S-matrix triggers bounds clearly.
        let s_active = vec![vec![c(0.2, 0.0), c(1.2, 0.0)], vec![c(0.4, 0.0), c(0.2, 0.0)]];
        let (sigma_max, recip_err) = s_matrix_passivity_and_reciprocity(&s_active);
        assert!(sigma_max > 1.0, "expected non-passive sigma_max > 1, got {sigma_max}");
        assert!(recip_err > 0.1, "expected non-reciprocal error > 0.1, got {recip_err}");
    }

    #[test]
    fn sparse_driven_solver_solve_many_checked_reports_backward_residual() {
        let a = CCsr::from_triplets(
            3,
            3,
            &[(0, 0, c(4.0, 0.5)), (1, 1, c(5.0, -0.2)), (2, 2, c(6.0, 0.1)), (0, 1, c(-1.0, 0.0)), (1, 0, c(-1.0, 0.0))],
        )
        .unwrap();
        let mut solver = SparseDrivenSolver::new("auto", 0.0);
        solver.update_matrix(&a).unwrap();
        let rhs = vec![vec![c(1.0, 0.0), c(0.0, 1.0), c(2.0, -1.0)]];
        let (xs, chk) = solver.solve_many_checked(&a, &rhs).unwrap();
        assert_eq!(xs.len(), 1);
        assert!(chk.algebraic_residual < 1e-13, "sparse driven residual = {:e}", chk.algebraic_residual);
    }
}
