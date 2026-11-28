// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Smoothed-aggregation algebraic multigrid for real and **complex symmetric** matrices
//! (port items "AMG for complex symmetric matrices", the multigrid half of the shifted-Laplacian
//! preconditioner P2, and FR-43a parallel Chebyshev & mixed-precision `AmgF32` cycles).
//!
//! * Strength of connection uses moduli relative to the largest off-diagonal of the row:
//!   `|a_ij| ≥ θ·max_k |a_ik|`.
//! * Aggregation: the standard three-pass greedy scheme; one real near-null-space vector
//!   (the constant) is carried down the hierarchy; the tentative prolongator is normalised per
//!   aggregate and then smoothed by one damped Jacobi step with the *complex* diagonal.
//! * Coarse operators are Galerkin products `R A P` with `R = Pᵀ` (plain transpose), which keeps a
//!   complex symmetric matrix complex symmetric. Set [`AmgOptions::hermitian`] to use `R = Pᴴ`.
//! * The smoother defaults to symmetric Gauss–Seidel (`SmootherKind::GaussSeidel`, forward before
//!   and backward after the coarse correction, preserving 0-ULP CPU goldens) or deterministic
//!   parallel Chebyshev polynomial smoothing (`SmootherKind::Chebyshev { degree }`) using a
//!   power-iteration spectral bound `λ_max(D⁻¹ A)` for multi-thread CPU and GPU parity.
//! * [`AmgF32`] runs the entire V-cycle hierarchy in single precision (`f32`) while implementing
//!   [`Preconditioner<f64>`] to drive outer double-precision PCG / FGMRES solves.

use crate::linalg::csr::Csr;
use crate::linalg::device::{CpuVec, DevicePreconditioner};
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::linop::Preconditioner;
use crate::linalg::par;
use crate::linalg::rng::random_vector;
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops::{self, CHUNK};
use std::cell::RefCell;

/// Multigrid level smoother selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SmootherKind {
    /// Sequential symmetric Gauss–Seidel (forward pre-smooth, backward post-smooth).
    #[default]
    GaussSeidel,
    /// Parallel deterministic Chebyshev polynomial smoother on `D⁻¹ A` of the given `degree`
    /// (typically `2..=4`), built purely from SpMVs and diagonal scalings.
    Chebyshev { degree: usize },
}

#[derive(Debug, Clone)]
pub struct AmgOptions {
    /// Strength threshold `θ` (relative to the row's largest off-diagonal modulus).
    pub theta: f64,
    /// Stop coarsening below this many unknowns.
    pub max_coarse: usize,
    pub max_levels: usize,
    /// Smooth the tentative prolongator (smoothed aggregation); `false` = plain aggregation.
    pub smooth_prolongator: bool,
    /// Gauss–Seidel sweeps before and after the coarse correction (when `smoother == GaussSeidel`).
    pub sweeps: usize,
    /// Use `R = Pᴴ` (Hermitian problems) instead of `R = Pᵀ` (complex symmetric).
    pub hermitian: bool,
    /// Smoother applied on each non-coarsest level (default [`SmootherKind::GaussSeidel`]).
    pub smoother: SmootherKind,
}

impl Default for AmgOptions {
    fn default() -> Self {
        Self {
            theta: 0.25,
            max_coarse: 150,
            max_levels: 20,
            smooth_prolongator: true,
            sweeps: 1,
            hermitian: false,
            smoother: SmootherKind::GaussSeidel,
        }
    }
}

struct Level<T: Scalar> {
    a: Csr<T>,
    inv_diag: Vec<T>,
    lmax: f64,
    /// Prolongation `P` (fine × coarse) and restriction `R` (coarse × fine); empty on the last level.
    p: Option<Csr<T>>,
    r: Option<Csr<T>>,
}

struct Scratch<T> {
    rhs: Vec<T>,
    x: Vec<T>,
    res: Vec<T>,
    tmp: Vec<T>,
}

/// Summary of the hierarchy, for logs and tests.
#[derive(Debug, Clone)]
pub struct AmgStats {
    pub level_sizes: Vec<usize>,
    pub level_nnz: Vec<usize>,
    /// `Σ nnz(A_l) / nnz(A_0)`.
    pub operator_complexity: f64,
}

pub struct Amg<T: FaerScalar> {
    levels: Vec<Level<T>>,
    coarse: ShiftedFactorization<T>,
    scratch: RefCell<Vec<Scratch<T>>>,
    opts: AmgOptions,
}

impl<T: FaerScalar> Amg<T> {
    pub fn new(a: &Csr<T>, opts: AmgOptions) -> Result<Self, SolverError> {
        if a.nrows() != a.ncols() {
            return Err(SolverError::InvalidArgument("AMG matrix must be square"));
        }
        let need_lmax = matches!(opts.smoother, SmootherKind::Chebyshev { .. });
        let raw_levels = build_hierarchy(a, &opts, need_lmax)?;
        let last = &raw_levels[raw_levels.len() - 1].a;
        let mut coarse = ShiftedFactorization::new(last, None, Definiteness::General)?;
        coarse.set_shift(T::ZERO)?;
        let scratch = raw_levels
            .iter()
            .map(|l| {
                let n = l.a.nrows();
                Scratch { rhs: vec![T::ZERO; n], x: vec![T::ZERO; n], res: vec![T::ZERO; n], tmp: vec![T::ZERO; n] }
            })
            .collect();
        Ok(Self { levels: raw_levels, coarse, scratch: RefCell::new(scratch), opts })
    }

    pub fn stats(&self) -> AmgStats {
        let level_sizes: Vec<usize> = self.levels.iter().map(|l| l.a.nrows()).collect();
        let level_nnz: Vec<usize> = self.levels.iter().map(|l| l.a.nnz()).collect();
        let denom = level_nnz.first().copied().unwrap_or(1).max(1) as f64;
        let operator_complexity = level_nnz.iter().sum::<usize>() as f64 / denom;
        AmgStats { level_sizes, level_nnz, operator_complexity }
    }

    fn cycle(&self, lvl: usize, sc: &mut [Scratch<T>]) {
        let Some((cur, rest)) = sc.split_first_mut() else {
            return;
        };
        let level = &self.levels[lvl];
        let (Some(p), Some(r_op)) = (&level.p, &level.r) else {
            cur.x.copy_from_slice(&cur.rhs);
            let _ = self.coarse.solve_block_in_place(&mut cur.x, 1);
            return;
        };
        cur.x.fill(T::ZERO);
        match self.opts.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.opts.sweeps {
                    gauss_seidel(&level.a, &level.inv_diag, &cur.rhs, &mut cur.x, true);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth(
                    &level.a,
                    &level.inv_diag,
                    level.lmax,
                    degree.max(1),
                    &cur.rhs,
                    &mut cur.x,
                    &mut cur.res,
                    &mut cur.tmp,
                    true,
                );
            }
        }
        // res = rhs − A x
        level.a.spmv(&cur.x, &mut cur.res);
        vecops::aypx(-T::ONE, &cur.rhs, &mut cur.res);
        r_op.spmv(&cur.res, &mut rest[0].rhs);
        self.cycle(lvl + 1, rest);
        p.spmv(&rest[0].x, &mut cur.tmp);
        vecops::axpy(T::ONE, &cur.tmp, &mut cur.x);
        match self.opts.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.opts.sweeps {
                    gauss_seidel(&level.a, &level.inv_diag, &cur.rhs, &mut cur.x, false);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth(
                    &level.a,
                    &level.inv_diag,
                    level.lmax,
                    degree.max(1),
                    &cur.rhs,
                    &mut cur.x,
                    &mut cur.res,
                    &mut cur.tmp,
                    false,
                );
            }
        }
    }
}

impl<T: FaerScalar> Preconditioner<T> for Amg<T> {
    /// One V-cycle from a zero initial guess: `z ≈ A⁻¹ r`.
    fn apply(&self, r: &[T], z: &mut [T]) {
        let mut sc = self.scratch.borrow_mut();
        sc[0].rhs.copy_from_slice(r);
        self.cycle(0, &mut sc);
        z.copy_from_slice(&sc[0].x);
    }
}

impl<T: FaerScalar> DevicePreconditioner<CpuVec<T>> for Amg<T> {
    #[inline]
    fn apply(&self, r: &CpuVec<T>, z: &mut CpuVec<T>) {
        Preconditioner::apply(self, &r.0, &mut z.0);
    }
}

fn build_hierarchy<T: Scalar>(a: &Csr<T>, opts: &AmgOptions, compute_lmax: bool) -> Result<Vec<Level<T>>, SolverError> {
    if a.values().iter().any(|v| !v.is_finite()) {
        return Err(SolverError::InvalidArgument("AMG matrix contains non-finite entries"));
    }
    let mut levels: Vec<Level<T>> = Vec::new();
    let mut cur = a.clone();
    let mut nullvec = vec![1.0f64; a.nrows()];
    loop {
        let n = cur.nrows();
        let inv_diag: Vec<T> =
            cur.diagonal().into_iter().map(|d| if d == T::ZERO { T::ONE } else { T::ONE / d }).collect();
        if n <= opts.max_coarse || levels.len() + 1 >= opts.max_levels {
            levels.push(Level { a: cur, inv_diag, lmax: 1.0, p: None, r: None });
            break;
        }
        let agg = aggregate(&cur, opts.theta);
        let nagg = agg.iter().copied().max().map_or(0, |m| m + 1);
        if nagg == 0 || nagg >= n {
            levels.push(Level { a: cur, inv_diag, lmax: 1.0, p: None, r: None });
            break;
        }
        let lmax = if compute_lmax { estimate_dinv_a_spectral_radius(&cur, &inv_diag) } else { 1.0 };
        let (p_tent, coarse_null) = tentative_prolongator(&agg, nagg, &nullvec)?;
        let p = if opts.smooth_prolongator { smooth(&cur, &inv_diag, &p_tent)? } else { p_tent };
        let r = {
            let t = p.transpose();
            if opts.hermitian { t.map(|v| v.conj()) } else { t }
        };
        let coarse = r.matmul(&cur.matmul(&p));
        levels.push(Level { a: cur, inv_diag, lmax, p: Some(p), r: Some(r) });
        cur = coarse;
        nullvec = coarse_null;
    }
    Ok(levels)
}

/// Deterministic power-iteration estimate of `λ_max(D⁻¹ A)` using the symmetrically scaled
/// operator `D⁻¹/² A D⁻¹/²` (20 steps with fixed seed and `1.05×` safety margin, falling back
/// to the Gershgorin row-sum bound if the Rayleigh vector vanishes).
pub fn estimate_dinv_a_spectral_radius<T: Scalar>(a: &Csr<T>, inv_diag: &[T]) -> f64 {
    let n = a.nrows();
    if n == 0 {
        return 1.0;
    }
    let gershgorin = || {
        (0..n)
            .map(|i| {
                let (_, v) = a.row(i);
                v.iter().map(|x| x.modulus()).sum::<f64>() * inv_diag[i].modulus()
            })
            .fold(0.0f64, f64::max)
            .max(1e-12)
    };
    let inv_sqrt: Vec<T> = inv_diag
        .iter()
        .map(|&d| {
            let s = d.sqrt_principal();
            if s.is_finite() { s } else { T::from_f64(d.modulus().sqrt()) }
        })
        .collect();
    let mut v = random_vector::<T>(n, 0xC4EB_5EED);
    let n0 = vecops::norm2(&v);
    if n0 <= 1e-300 || !n0.is_finite() {
        return gershgorin();
    }
    vecops::scale_real(1.0 / n0, &mut v);
    let mut sv = vec![T::ZERO; n];
    let mut w = vec![T::ZERO; n];
    let mut est = 0.0f64;
    let par = par::worth_parallel(n);
    for _ in 0..20 {
        par::for_chunks_mut(&mut sv, CHUNK, par, |c, svc| {
            let s = c * CHUNK;
            let len = svc.len();
            for (( dst, &vi), &si) in svc.iter_mut().zip(&v[s..s + len]).zip(&inv_sqrt[s..s + len]) {
                *dst = si * vi;
            }
        });
        a.spmv(&sv, &mut w);
        par::for_chunks_mut(&mut w, CHUNK, par, |c, wc| {
            let s = c * CHUNK;
            let len = wc.len();
            for (wi, &si) in wc.iter_mut().zip(&inv_sqrt[s..s + len]) {
                *wi *= si;
            }
        });
        let nrm = vecops::norm2(&w);
        if nrm <= 1e-300 || !nrm.is_finite() {
            return gershgorin();
        }
        est = nrm;
        vecops::scale_real(1.0 / nrm, &mut w);
        std::mem::swap(&mut v, &mut w);
    }
    (est * 1.05).max(1e-12)
}

/// Deterministic parallel Chebyshev polynomial smoother on `D⁻¹ A` over `[1.1·lmax / 30, 1.1·lmax]`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn chebyshev_smooth<T: Scalar>(
    a: &Csr<T>,
    inv_diag: &[T],
    lmax: f64,
    degree: usize,
    b: &[T],
    x: &mut [T],
    res: &mut [T],
    d: &mut [T],
    zero_guess: bool,
) {
    let n = x.len();
    if n == 0 || degree == 0 {
        return;
    }
    let hi = 1.1 * lmax.max(1e-12);
    let lo = hi / 30.0;
    let theta = 0.5 * (hi + lo);
    let delta = 0.5 * (hi - lo);
    let sigma = theta / delta;
    let inv_theta = 1.0 / theta;
    let par = par::worth_parallel(n);

    if zero_guess {
        par::for_chunks_mut(d, CHUNK, par, |c, dc| {
            let s = c * CHUNK;
            for (i, di) in dc.iter_mut().enumerate() {
                *di = (inv_diag[s + i] * b[s + i]).scale(inv_theta);
            }
        });
        x.copy_from_slice(d);
    } else {
        a.spmv(x, res);
        vecops::aypx(-T::ONE, b, res);
        par::for_chunks_mut(d, CHUNK, par, |c, dc| {
            let s = c * CHUNK;
            for (i, di) in dc.iter_mut().enumerate() {
                *di = (inv_diag[s + i] * res[s + i]).scale(inv_theta);
            }
        });
        vecops::axpy(T::ONE, d, x);
    }

    let mut rho_old = 1.0 / sigma;
    for _ in 1..degree {
        let rho_new = 1.0 / (2.0 * sigma - rho_old);
        let c1 = rho_new * rho_old;
        let c2 = 2.0 * rho_new / delta;
        a.spmv(x, res);
        vecops::aypx(-T::ONE, b, res);
        par::for_chunks_mut(d, CHUNK, par, |c, dc| {
            let s = c * CHUNK;
            for (i, di) in dc.iter_mut().enumerate() {
                *di = di.scale(c1) + (inv_diag[s + i] * res[s + i]).scale(c2);
            }
        });
        vecops::axpy(T::ONE, d, x);
        rho_old = rho_new;
    }
}

/// One Gauss–Seidel sweep (`forward` or backward) for `A x = b`, in place on `x`.
pub(crate) fn gauss_seidel<T: Scalar>(a: &Csr<T>, inv_diag: &[T], b: &[T], x: &mut [T], forward: bool) {
    let n = a.nrows();
    let mut sweep = |i: usize| {
        let (cols, vals) = a.row(i);
        let mut s = b[i];
        for (&j, &v) in cols.iter().zip(vals) {
            s -= v * x[j as usize];
        }
        // s = b_i − Σ_j a_ij x_j (diagonal included); add back the diagonal contribution.
        x[i] += s * inv_diag[i];
    };
    if forward {
        (0..n).for_each(&mut sweep);
    } else {
        (0..n).rev().for_each(&mut sweep);
    }
}

// ---------------------------------------------------------------------------------------
// Mixed-Precision Single-Precision V-Cycle Hierarchy (`AmgF32`)
// ---------------------------------------------------------------------------------------

/// Compact single-precision `f32` CSR matrix used by [`AmgF32`] and `AmsF32`.
#[derive(Debug, Clone)]
pub(crate) struct CsrF32 {
    pub(crate) nrows: usize,
    pub(crate) ncols: usize,
    pub(crate) indptr: Vec<usize>,
    pub(crate) indices: Vec<u32>,
    pub(crate) values: Vec<f32>,
}

impl CsrF32 {
    pub(crate) fn from_csr_f64(a: &Csr<f64>) -> Self {
        Self {
            nrows: a.nrows(),
            ncols: a.ncols(),
            indptr: a.indptr().to_vec(),
            indices: a.indices().to_vec(),
            values: a.values().iter().map(|&v| v as f32).collect(),
        }
    }

    #[inline]
    pub(crate) fn nnz(&self) -> usize {
        self.values.len()
    }

    #[inline]
    pub(crate) fn row(&self, i: usize) -> (&[u32], &[f32]) {
        let r = self.indptr[i]..self.indptr[i + 1];
        (&self.indices[r.clone()], &self.values[r])
    }

    pub(crate) fn spmv(&self, x: &[f32], y: &mut [f32]) {
        if x.len() != self.ncols || y.len() != self.nrows {
            return;
        }
        let par = par::worth_parallel(self.values.len());
        par::for_chunks_mut(y, CHUNK, par, |c, yc| {
            let start = c * CHUNK;
            let ip = &self.indptr[start..=start + yc.len()];
            for (yi, w) in yc.iter_mut().zip(ip.windows(2)) {
                let mut s = 0.0f32;
                for (&col, &val) in self.indices[w[0]..w[1]].iter().zip(&self.values[w[0]..w[1]]) {
                    s += val * x[col as usize];
                }
                *yi = s;
            }
        });
    }
}

pub(crate) fn gauss_seidel_f32(a: &CsrF32, inv_diag: &[f32], b: &[f32], x: &mut [f32], forward: bool) {
    let n = a.nrows;
    let mut sweep = |i: usize| {
        let (cols, vals) = a.row(i);
        let mut s = b[i];
        for (&j, &v) in cols.iter().zip(vals) {
            s -= v * x[j as usize];
        }
        x[i] += s * inv_diag[i];
    };
    if forward {
        (0..n).for_each(&mut sweep);
    } else {
        (0..n).rev().for_each(&mut sweep);
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn chebyshev_smooth_f32(
    a: &CsrF32,
    inv_diag: &[f32],
    lmax: f32,
    degree: usize,
    b: &[f32],
    x: &mut [f32],
    res: &mut [f32],
    d: &mut [f32],
    zero_guess: bool,
) {
    let n = x.len();
    if n == 0 || degree == 0 {
        return;
    }
    let hi = 1.1f32 * lmax.max(1e-12);
    let lo = hi / 30.0f32;
    let theta = 0.5f32 * (hi + lo);
    let delta = 0.5f32 * (hi - lo);
    let sigma = theta / delta;
    let inv_theta = 1.0f32 / theta;
    let par = par::worth_parallel(n);

    if zero_guess {
        par::for_chunks_mut(d, CHUNK, par, |c, dc| {
            let s = c * CHUNK;
            for (i, di) in dc.iter_mut().enumerate() {
                *di = inv_diag[s + i] * b[s + i] * inv_theta;
            }
        });
        x.copy_from_slice(d);
    } else {
        a.spmv(x, res);
        par::for_chunks_mut(res, CHUNK, par, |c, rc| {
            let s = c * CHUNK;
            for (i, ri) in rc.iter_mut().enumerate() {
                *ri = b[s + i] - *ri;
            }
        });
        par::for_chunks_mut(d, CHUNK, par, |c, dc| {
            let s = c * CHUNK;
            for (i, di) in dc.iter_mut().enumerate() {
                *di = inv_diag[s + i] * res[s + i] * inv_theta;
            }
        });
        par::for_chunks_mut(x, CHUNK, par, |c, xc| {
            let s = c * CHUNK;
            for (i, xi) in xc.iter_mut().enumerate() {
                *xi += d[s + i];
            }
        });
    }

    let mut rho_old = 1.0f32 / sigma;
    for _ in 1..degree {
        let rho_new = 1.0f32 / (2.0f32 * sigma - rho_old);
        let c1 = rho_new * rho_old;
        let c2 = 2.0f32 * rho_new / delta;
        a.spmv(x, res);
        par::for_chunks_mut(res, CHUNK, par, |c, rc| {
            let s = c * CHUNK;
            for (i, ri) in rc.iter_mut().enumerate() {
                *ri = b[s + i] - *ri;
            }
        });
        par::for_chunks_mut(d, CHUNK, par, |c, dc| {
            let s = c * CHUNK;
            for (i, di) in dc.iter_mut().enumerate() {
                *di = c1 * *di + c2 * inv_diag[s + i] * res[s + i];
            }
        });
        par::for_chunks_mut(x, CHUNK, par, |c, xc| {
            let s = c * CHUNK;
            for (i, xi) in xc.iter_mut().enumerate() {
                *xi += d[s + i];
            }
        });
        rho_old = rho_new;
    }
}

pub(crate) struct LevelF32 {
    #[cfg(feature = "gpu")]
    pub(crate) a_f64: Csr<f64>,
    pub(crate) a: CsrF32,
    pub(crate) inv_diag: Vec<f32>,
    pub(crate) lmax: f32,
    #[cfg(feature = "gpu")]
    pub(crate) p_f64: Option<Csr<f64>>,
    #[cfg(feature = "gpu")]
    pub(crate) r_f64: Option<Csr<f64>>,
    pub(crate) p: Option<CsrF32>,
    pub(crate) r: Option<CsrF32>,
}

struct ScratchF32 {
    rhs: Vec<f32>,
    x: Vec<f32>,
    res: Vec<f32>,
    tmp: Vec<f32>,
}

/// Mixed-precision smoothed-aggregation AMG preconditioner: stores and executes the V-cycle
/// hierarchy in `f32` while exposing [`Preconditioner<f64>`] for outer `f64` PCG / FGMRES solves.
pub struct AmgF32 {
    levels: Vec<LevelF32>,
    coarse: ShiftedFactorization<f64>,
    coarse_buf: RefCell<Vec<f64>>,
    scratch: RefCell<Vec<ScratchF32>>,
    opts: AmgOptions,
}

impl AmgF32 {
    pub fn new(a: &Csr<f64>, opts: AmgOptions) -> Result<Self, SolverError> {
        if a.nrows() != a.ncols() {
            return Err(SolverError::InvalidArgument("AmgF32 matrix must be square"));
        }
        // Always compute spectral bounds so GPU / Chebyshev callers have `lmax` ready.
        let raw = build_hierarchy(a, &opts, true)?;
        let last_a = &raw[raw.len() - 1].a;
        let mut coarse = ShiftedFactorization::new(last_a, None, Definiteness::General)?;
        coarse.set_shift(0.0)?;
        let n_coarse = last_a.nrows();

        let levels: Vec<LevelF32> = raw
            .into_iter()
            .map(|l| LevelF32 {
                a: CsrF32::from_csr_f64(&l.a),
                #[cfg(feature = "gpu")]
                a_f64: l.a,
                inv_diag: l.inv_diag.iter().map(|&d| d as f32).collect(),
                lmax: l.lmax as f32,
                p: l.p.as_ref().map(CsrF32::from_csr_f64),
                r: l.r.as_ref().map(CsrF32::from_csr_f64),
                #[cfg(feature = "gpu")]
                p_f64: l.p,
                #[cfg(feature = "gpu")]
                r_f64: l.r,
            })
            .collect();

        let scratch = levels
            .iter()
            .map(|l| {
                let n = l.a.nrows;
                ScratchF32 { rhs: vec![0.0; n], x: vec![0.0; n], res: vec![0.0; n], tmp: vec![0.0; n] }
            })
            .collect();

        Ok(Self {
            levels,
            coarse,
            coarse_buf: RefCell::new(vec![0.0; n_coarse]),
            scratch: RefCell::new(scratch),
            opts,
        })
    }

    /// Convenience constructor building `AmgF32` with a deterministic parallel Chebyshev smoother
    /// of the given `degree` (e.g. `3`).
    pub fn new_chebyshev(a: &Csr<f64>, degree: usize) -> Result<Self, SolverError> {
        Self::new(
            a,
            AmgOptions { smoother: SmootherKind::Chebyshev { degree: degree.max(1) }, ..Default::default() },
        )
    }

    pub fn dim(&self) -> usize {
        self.levels.first().map_or(0, |l| l.a.nrows)
    }

    pub fn stats(&self) -> AmgStats {
        let level_sizes: Vec<usize> = self.levels.iter().map(|l| l.a.nrows).collect();
        let level_nnz: Vec<usize> = self.levels.iter().map(|l| l.a.nnz()).collect();
        let denom = level_nnz.first().copied().unwrap_or(1).max(1) as f64;
        let operator_complexity = level_nnz.iter().sum::<usize>() as f64 / denom;
        AmgStats { level_sizes, level_nnz, operator_complexity }
    }

    #[cfg(feature = "gpu")]
    pub(crate) fn levels_f32(&self) -> &[LevelF32] {
        &self.levels
    }

    pub fn options(&self) -> &AmgOptions {
        &self.opts
    }

    pub(crate) fn solve_coarse_f32(&self, rhs: &[f32], x: &mut [f32]) {
        let mut buf = self.coarse_buf.borrow_mut();
        for (dst, &src) in buf.iter_mut().zip(rhs) {
            *dst = src as f64;
        }
        let _ = self.coarse.solve_block_in_place(&mut buf, 1);
        for (dst, &src) in x.iter_mut().zip(buf.iter()) {
            *dst = src as f32;
        }
    }

    /// Runs one `f32` V-cycle directly on `f32` slices.
    pub fn apply_f32(&self, r: &[f32], z: &mut [f32]) {
        let mut sc = self.scratch.borrow_mut();
        sc[0].rhs.copy_from_slice(r);
        self.cycle_f32(0, &mut sc);
        z.copy_from_slice(&sc[0].x);
    }

    fn cycle_f32(&self, lvl: usize, sc: &mut [ScratchF32]) {
        let Some((cur, rest)) = sc.split_first_mut() else {
            return;
        };
        let level = &self.levels[lvl];
        let (Some(p), Some(r_op)) = (&level.p, &level.r) else {
            self.solve_coarse_f32(&cur.rhs, &mut cur.x);
            return;
        };
        let n = level.a.nrows;
        let par = par::worth_parallel(n);
        cur.x.fill(0.0);
        match self.opts.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.opts.sweeps {
                    gauss_seidel_f32(&level.a, &level.inv_diag, &cur.rhs, &mut cur.x, true);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth_f32(
                    &level.a,
                    &level.inv_diag,
                    level.lmax,
                    degree.max(1),
                    &cur.rhs,
                    &mut cur.x,
                    &mut cur.res,
                    &mut cur.tmp,
                    true,
                );
            }
        }
        level.a.spmv(&cur.x, &mut cur.res);
        par::for_chunks_mut(&mut cur.res, CHUNK, par, |c, rc| {
            let s = c * CHUNK;
            for (i, ri) in rc.iter_mut().enumerate() {
                *ri = cur.rhs[s + i] - *ri;
            }
        });
        r_op.spmv(&cur.res, &mut rest[0].rhs);
        self.cycle_f32(lvl + 1, rest);
        p.spmv(&rest[0].x, &mut cur.tmp);
        par::for_chunks_mut(&mut cur.x, CHUNK, par, |c, xc| {
            let s = c * CHUNK;
            for (i, xi) in xc.iter_mut().enumerate() {
                *xi += cur.tmp[s + i];
            }
        });
        match self.opts.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.opts.sweeps {
                    gauss_seidel_f32(&level.a, &level.inv_diag, &cur.rhs, &mut cur.x, false);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth_f32(
                    &level.a,
                    &level.inv_diag,
                    level.lmax,
                    degree.max(1),
                    &cur.rhs,
                    &mut cur.x,
                    &mut cur.res,
                    &mut cur.tmp,
                    false,
                );
            }
        }
    }
}

/// Computes a power-of-two scaling `(to_f32_scale, from_f32_scale)` so extreme `f64` residuals
/// (`‖r‖ < 1e-20` or `‖r‖ > 1e20`) never underflow or overflow in `f32`, while normal residuals
/// (`1e-20 <= ‖r‖ <= 1e20`) use exact `1.0` scaling with zero extra rounding.
#[inline]
pub(crate) fn f32_safe_scale(r: &[f64]) -> (f64, f64) {
    let nrm = vecops::norm2(r);
    if nrm == 0.0 || !nrm.is_finite() || (1e-20..=1e20).contains(&nrm) {
        (1.0, 1.0)
    } else {
        let exp = nrm.log2().floor() as i32;
        let inv = 2.0f64.powi(-exp);
        let fwd = 2.0f64.powi(exp);
        (inv, fwd)
    }
}

impl Preconditioner<f64> for AmgF32 {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        let (to_s, from_s) = f32_safe_scale(r);
        let mut sc = self.scratch.borrow_mut();
        for (dst, &src) in sc[0].rhs.iter_mut().zip(r) {
            *dst = (src * to_s) as f32;
        }
        self.cycle_f32(0, &mut sc);
        for (dst, &src) in z.iter_mut().zip(&sc[0].x) {
            *dst = (src as f64) * from_s;
        }
    }
}

impl DevicePreconditioner<CpuVec<f64>> for AmgF32 {
    #[inline]
    fn apply(&self, r: &CpuVec<f64>, z: &mut CpuVec<f64>) {
        Preconditioner::apply(self, &r.0, &mut z.0);
    }
}

/// Standard three-pass aggregation on the strength graph; returns the aggregate of each node.
fn aggregate<T: Scalar>(a: &Csr<T>, theta: f64) -> Vec<usize> {
    let n = a.nrows();
    let rowmax: Vec<f64> = (0..n)
        .map(|i| {
            let (cols, vals) = a.row(i);
            cols.iter().zip(vals).filter(|&(&j, _)| j as usize != i).map(|(_, v)| v.modulus()).fold(0.0, f64::max)
        })
        .collect();
    let strong = |i: usize| -> Vec<(usize, f64)> {
        let (cols, vals) = a.row(i);
        cols.iter()
            .zip(vals)
            .filter(|&(&j, v)| j as usize != i && v.modulus() > 0.0 && v.modulus() >= theta * rowmax[i])
            .map(|(&j, v)| (j as usize, v.modulus()))
            .collect()
    };
    let nbrs: Vec<Vec<(usize, f64)>> = (0..n).map(strong).collect();
    const NONE: usize = usize::MAX;
    let mut agg = vec![NONE; n];
    let mut next = 0usize;
    let mut blocked = vec![false; n];
    for i in 0..n {
        if !blocked[i] && !nbrs[i].is_empty() {
            blocked[i] = true;
            agg[i] = next;
            for &(j, _) in &nbrs[i] {
                blocked[j] = true;
                if agg[j] == NONE {
                    agg[j] = next;
                }
                for &(k, _) in &nbrs[j] {
                    blocked[k] = true;
                }
            }
            next += 1;
        }
    }
    loop {
        let snapshot = agg.clone();
        let mut changed = false;
        for i in 0..n {
            if agg[i] != NONE {
                continue;
            }
            if let Some(&(j, _)) =
                nbrs[i].iter().filter(|&&(j, _)| snapshot[j] != NONE).max_by(|x, y| x.1.total_cmp(&y.1))
            {
                agg[i] = snapshot[j];
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    for a in agg.iter_mut() {
        if *a == NONE {
            *a = next;
            next += 1;
        }
    }
    agg
}

/// `P_tent` with columns normalised per aggregate for the near-null-space vector `b`; also
/// returns the coarse near-null-space vector (the aggregate norms).
fn tentative_prolongator<T: Scalar>(agg: &[usize], nagg: usize, b: &[f64]) -> Result<(Csr<T>, Vec<f64>), SolverError> {
    let mut nrm = vec![0.0f64; nagg];
    for (i, &g) in agg.iter().enumerate() {
        nrm[g] += b[i] * b[i];
    }
    for v in &mut nrm {
        *v = v.sqrt();
    }
    let n = agg.len();
    let indptr: Vec<usize> = (0..=n).collect();
    let indices: Vec<u32> = agg.iter().map(|&g| g as u32).collect();
    let values: Vec<T> =
        (0..n).map(|i| T::from_f64(if nrm[agg[i]] > 0.0 { b[i] / nrm[agg[i]] } else { 0.0 })).collect();
    let p = Csr::try_new(n, nagg, indptr, indices, values)
        .map_err(|_| SolverError::InvalidArgument("failed to construct tentative prolongator"))?;
    Ok((p, nrm))
}

/// `P = (I − ω D⁻¹ A) P_tent` with `ω = 4/(3ρ)` and `ρ` the Gershgorin bound of `D⁻¹A` in moduli.
fn smooth<T: Scalar>(a: &Csr<T>, inv_diag: &[T], p_tent: &Csr<T>) -> Result<Csr<T>, SolverError> {
    let n = a.nrows();
    let rho = (0..n)
        .map(|i| {
            let (_, v) = a.row(i);
            v.iter().map(|x| x.modulus()).sum::<f64>() * inv_diag[i].modulus()
        })
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);
    let omega = 4.0 / (3.0 * rho);
    let ap = a.matmul(p_tent);
    let scaled = Csr::try_new(
        ap.nrows(),
        ap.ncols(),
        ap.indptr().to_vec(),
        ap.indices().to_vec(),
        (0..n)
            .flat_map(|i| {
                let (_, v) = ap.row(i);
                v.iter().map(move |x| *x * inv_diag[i].scale(-omega)).collect::<Vec<_>>()
            })
            .collect(),
    )
    .map_err(|_| SolverError::InvalidArgument("failed to construct smoothed prolongator"))?;
    Ok(Csr::axpby(T::ONE, p_tent, T::ONE, &scaled))
}

/// Builds the shifted-Laplacian operator `K − (β₁ + iβ₂)·k²·M` (Erlangga et al.): the matrix that
/// is *inverted approximately* (one AMG V-cycle) to precondition the Helmholtz/PML operator
/// `K − k²M (+ damping)`. The usual choice is `β = 1 + 0.5i`; larger imaginary parts make the
/// multigrid easier and the preconditioner a worse approximation.
pub fn shifted_laplacian_operator<T: Scalar>(
    k: &Csr<T>,
    mass: Option<&Csr<T>>,
    k2: f64,
    beta: crate::linalg::scalar::C64,
) -> Csr<crate::linalg::scalar::C64> {
    let kc = k.to_complex();
    let shift = beta.scale(-k2);
    match mass {
        Some(m) => Csr::axpby(crate::linalg::scalar::C64::ONE, &kc, shift, &m.to_complex()),
        None => Csr::axpby(crate::linalg::scalar::C64::ONE, &kc, shift, &Csr::identity(k.nrows())),
    }
}
