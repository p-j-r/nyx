//! Chebyshev-filtered subspace iteration for many lowest eigenpairs (port item E5; the idea of ChASE,
//! Winkelmann et al. 2019): for `nev` in the tens to hundreds (photonic-crystal band structures)
//! where a polynomial filter beats Krylov methods because it needs only matrix products and block
//! orthogonalisation, and the work parallelises over the block.
//!
//! Each outer iteration applies a Chebyshev polynomial of degree `m` that damps the unwanted interval `[c, b]`
//! (`b` an upper bound of the spectrum from Lanczos, `c` the cutoff: the largest Ritz value of the block,
//! so the wanted eigenvalues lie below it) with the Zhou–Saad scaling that keeps the polynomial values at the lower
//! end of the spectrum `O(1)`; then orthonormalises, does Rayleigh–Ritz and locks converged vectors (they
//! are no longer filtered). Standard Hermitian problems only (`B = I`); a generalised problem needs the
//! Cholesky factor of `B` (not implemented).

use crate::solvers::csr::Csr;
use crate::solvers::dense::{Dense, DenseError, hermitian_eig};
use crate::solvers::faer_glue::FaerScalar;
use crate::solvers::rng::random_vector;
use crate::solvers::scalar::Scalar;
use crate::solvers::vecops::{dotc, norm2};

#[derive(Debug, Clone)]
pub struct ChebOptions {
    pub nev: usize,
    /// Extra vectors in the block (default `max(nev/4, 4)`): the cutoff is the largest Ritz value of the block.
    pub nex: usize,
    /// Polynomial degree per outer iteration.
    pub degree: usize,
    pub tol: f64,
    pub max_it: usize,
    /// Lanczos steps for the spectral upper bound.
    pub lanczos_steps: usize,
    pub seed: u64,
}

impl ChebOptions {
    pub fn new(nev: usize) -> Self {
        Self { nev, nex: (nev / 4).max(4), degree: 20, tol: 1e-8, max_it: 100, lanczos_steps: 30, seed: 0xC4EB }
    }
}

#[derive(Debug, Clone)]
pub struct ChebResult<T: Scalar> {
    pub values: Vec<f64>,
    pub vectors: Vec<Vec<T>>,
    /// Relative residuals `‖Ax − λx‖ / max(|λ|, ε‖A‖)`.
    pub residuals: Vec<f64>,
    pub iterations: usize,
    pub matvecs: usize,
    pub converged: bool,
    /// Number of locked vectors after each iteration.
    pub locked_history: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ChebError {
    InvalidArgument(&'static str),
    Dense(DenseError),
}

impl std::fmt::Display for ChebError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ChebError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
            ChebError::Dense(e) => write!(f, "dense step failed: {e}"),
        }
    }
}

impl std::error::Error for ChebError {}
impl From<DenseError> for ChebError {
    fn from(e: DenseError) -> Self {
        ChebError::Dense(e)
    }
}

/// Upper bound of the spectrum of the Hermitian `a`: largest Lanczos Ritz value plus its residual (a rigorous bound for Hermitian matrices).
fn spectral_upper_bound<T: FaerScalar>(a: &Csr<T>, steps: usize, seed: u64) -> Result<f64, DenseError> {
    let n = a.nrows();
    let k = steps.min(n - 1).max(2);
    let mut v = random_vector::<T>(n, seed);
    let nv = norm2(&v);
    v.iter_mut().for_each(|x| *x = x.scale(1.0 / nv));
    let mut vprev = vec![T::ZERO; n];
    let (mut alpha, mut beta) = (Vec::new(), Vec::new());
    let mut w = vec![T::ZERO; n];
    let mut b_prev = 0.0;
    let mut last_beta = 0.0;
    for _ in 0..k {
        a.spmv(&v, &mut w);
        let al = dotc(&v, &w).re();
        for i in 0..n {
            w[i] = w[i] - v[i].scale(al) - vprev[i].scale(b_prev);
        }
        alpha.push(al);
        let bt = norm2(&w);
        last_beta = bt;
        if bt < 1e-12 {
            break;
        }
        beta.push(bt);
        vprev = v.clone();
        v = w.iter().map(|x| x.scale(1.0 / bt)).collect();
        b_prev = bt;
    }
    let m = alpha.len();
    let t = Dense::<f64>::from_fn(m, m, |i, j| {
        if i == j {
            alpha[i]
        } else if i + 1 == j || j + 1 == i {
            beta[i.min(j)]
        } else {
            0.0
        }
    });
    let (w, q) = hermitian_eig(&t)?;
    // largest Ritz value + |β_m · last component of its eigenvector|
    let top = w[m - 1] + last_beta * q[(m - 1, m - 1)].abs();
    Ok(top)
}

fn orthonormalise<T: Scalar>(block: &mut Vec<Vec<T>>, locked: &[Vec<T>]) {
    let mut out: Vec<Vec<T>> = Vec::with_capacity(block.len());
    for mut v in block.drain(..) {
        for _ in 0..2 {
            for u in locked.iter().chain(out.iter()) {
                let h = dotc(u, &v);
                for (vi, ui) in v.iter_mut().zip(u) {
                    *vi -= h * *ui;
                }
            }
        }
        let nv = norm2(&v);
        if nv > 1e-12 {
            v.iter_mut().for_each(|x| *x = x.scale(1.0 / nv));
            out.push(v);
        }
    }
    *block = out;
}

pub fn chebyshev_subspace<T: FaerScalar>(a: &Csr<T>, opts: &ChebOptions) -> Result<ChebResult<T>, ChebError> {
    let n = a.nrows();
    let m = opts.nev + opts.nex;
    if opts.nev == 0 || a.ncols() != n || 2 * m >= n || opts.degree < 2 {
        return Err(ChebError::InvalidArgument("need nev ≥ 1, a square matrix, 2(nev + nex) < n and degree ≥ 2"));
    }
    let bsup = spectral_upper_bound(a, opts.lanczos_steps, opts.seed ^ 0x55)?;
    let mut matvecs = opts.lanczos_steps.min(n);
    // starting block and a first Rayleigh–Ritz (no filtering) for the initial cutoff and lower bound
    let mut block: Vec<Vec<T>> =
        (0..m).map(|j| random_vector::<T>(n, opts.seed.wrapping_add(j as u64 * 7919))).collect();
    orthonormalise(&mut block, &[]);
    let mut locked: Vec<Vec<T>> = Vec::new();
    let mut locked_vals: Vec<f64> = Vec::new();
    let mut locked_res: Vec<f64> = Vec::new();
    let mut locked_history = Vec::new();
    let mut theta: Vec<f64> = Vec::new();
    let mut ablock: Vec<Vec<T>> = Vec::new();
    let rayleigh_ritz = |block: &mut Vec<Vec<T>>,
                         ablock: &mut Vec<Vec<T>>,
                         theta: &mut Vec<f64>,
                         matvecs: &mut usize|
     -> Result<(), ChebError> {
        let k = block.len();
        *ablock = block
            .iter()
            .map(|v| {
                let mut y = vec![T::ZERO; n];
                a.spmv(v, &mut y);
                y
            })
            .collect();
        *matvecs += k;
        let g =
            Dense::from_fn(k, k, |i, j| (dotc(&block[i], &ablock[j]) + dotc(&ablock[i], &block[j]).conj()).scale(0.5));
        let (w, q) = hermitian_eig(&g)?;
        let rot = |vecs: &Vec<Vec<T>>| -> Vec<Vec<T>> {
            (0..k)
                .map(|c| {
                    let mut out = vec![T::ZERO; n];
                    for (i, v) in vecs.iter().enumerate() {
                        let coef = q[(i, c)];
                        for (o, vi) in out.iter_mut().zip(v) {
                            *o += coef * *vi;
                        }
                    }
                    out
                })
                .collect()
        };
        let nb = rot(block);
        let na = rot(ablock);
        *block = nb;
        *ablock = na;
        *theta = w;
        Ok(())
    };
    rayleigh_ritz(&mut block, &mut ablock, &mut theta, &mut matvecs)?;
    let mut iterations = 0;
    let mut converged = false;
    let anorm = bsup.abs().max(1e-300);
    for it in 0..opts.max_it {
        iterations = it + 1;
        // residuals of the active block; lock converged vectors from the bottom of the spectrum upwards
        let mut nlock_now = 0;
        for j in 0..block.len() {
            if locked.len() + j >= opts.nev && j > 0 {
                break;
            }
            let mut rj = ablock[j].clone();
            for (ri, vi) in rj.iter_mut().zip(&block[j]) {
                *ri -= vi.scale(theta[j]);
            }
            let res = norm2(&rj) / theta[j].abs().max(1e-8 * anorm);
            if res <= opts.tol {
                nlock_now += 1;
                locked_res.push(res);
            } else {
                locked_res.truncate(locked.len());
                break;
            }
        }
        for _ in 0..nlock_now {
            locked.push(block.remove(0));
            ablock.remove(0);
            locked_vals.push(theta.remove(0));
        }
        locked_history.push(locked.len());
        if locked.len() >= opts.nev {
            converged = true;
            break;
        }
        // Chebyshev filter on the active block: damp [c, b], c = largest Ritz value of the block
        let c = theta[theta.len() - 1];
        let a0 = theta[0].min(locked_vals.first().copied().unwrap_or(theta[0])) - 1e-3 * anorm;
        let (e, h) = (0.5 * (bsup + c), 0.5 * (bsup - c));
        if h <= 0.0 {
            break;
        }
        let sigma1 = h / (a0 - e);
        let apply_shifted = |x: &Vec<T>, out: &mut Vec<T>| {
            a.spmv(x, out);
            for (o, xi) in out.iter_mut().zip(x) {
                *o -= xi.scale(e);
            }
        };
        let mut filtered: Vec<Vec<T>> = Vec::with_capacity(block.len());
        for x in &block {
            // Y₁ = (σ₁/h)(A − eI) x
            let mut y = vec![T::ZERO; n];
            apply_shifted(x, &mut y);
            y.iter_mut().for_each(|v| *v = v.scale(sigma1 / h));
            matvecs += 1;
            let (mut xm, mut sigma) = (x.clone(), sigma1);
            let mut tmp = vec![T::ZERO; n];
            for _ in 2..=opts.degree {
                let sigma_new = 1.0 / (2.0 / sigma1 - sigma);
                apply_shifted(&y, &mut tmp);
                matvecs += 1;
                let ynew: Vec<T> = tmp
                    .iter()
                    .zip(&xm)
                    .map(|(t, xv)| t.scale(2.0 * sigma_new / h) - xv.scale(sigma * sigma_new))
                    .collect();
                xm = std::mem::replace(&mut y, ynew);
                sigma = sigma_new;
            }
            filtered.push(y);
        }
        block = filtered;
        orthonormalise(&mut block, &locked);
        if block.len() + locked.len() < opts.nev {
            return Err(ChebError::InvalidArgument("the filtered block lost rank: increase nex or reduce the degree"));
        }
        rayleigh_ritz(&mut block, &mut ablock, &mut theta, &mut matvecs)?;
    }
    // assemble: locked vectors, then the best of the active block up to nev
    let mut values = locked_vals.clone();
    let mut vectors = locked.clone();
    let mut residuals = locked_res.clone();
    residuals.resize(locked.len(), 0.0);
    for j in 0..block.len() {
        if values.len() >= opts.nev {
            break;
        }
        let mut rj = ablock[j].clone();
        for (ri, vi) in rj.iter_mut().zip(&block[j]) {
            *ri -= vi.scale(theta[j]);
        }
        residuals.push(norm2(&rj) / theta[j].abs().max(1e-8 * anorm));
        values.push(theta[j]);
        vectors.push(block[j].clone());
    }
    values.truncate(opts.nev);
    vectors.truncate(opts.nev);
    residuals.truncate(opts.nev);
    Ok(ChebResult { values, vectors, residuals, iterations, matvecs, converged, locked_history })
}
