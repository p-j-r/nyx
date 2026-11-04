//! The entry points the physics kernels call: a direct solve of `A x = b`, eigenpairs of `A x = lambda B x` near a shift, and small
//! dense eigenproblems. Everything else of `solvers` (Krylov methods, preconditioners, factorisation variants) stays reachable
//! through its own module; the kernels use only this file.
//!
//! What the driver adds to the library:
//! * **progress and cancel** (`crate::progress`): a stage line per factorisation and per eigensolve, and a cancel flag checked at
//!   every operator application, so a window can stop a long eigensolve (the library itself knows nothing about `Progress`);
//! * the **factorisation pool** (`crate::par::on_factor_pool`): faer factorisations scale to about 8 threads only;
//! * **projection** of a known null space out of every operator application (the discrete gradients of the curl-curl operator);
//! * the **real non-symmetric** pencil of the 2D waveguide modes, which the library does not take as real (it promotes to complex
//!   and keeps the real eigenvalues).
//!
//! All errors are plain-language strings, as everywhere in the kernels.

use crate::progress;
use crate::solvers::csr::Csr;
use crate::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::solvers::eps::{EpsOptions, EpsReason, Kind, Which, krylov_schur};
use crate::solvers::faer_glue::FaerScalar;
use crate::solvers::linop::LinearOperator;
use crate::solvers::problem::Symmetry;
use crate::solvers::scalar::{C64, Scalar};
use crate::solvers::vecops::norm2;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::panic::{AssertUnwindSafe, catch_unwind, resume_unwind};

// ------------------------------------------------------------------------------------------------------------ direct solves

/// A sparse matrix factorised once, then solved for any number of right-hand sides: Cholesky when the matrix is symmetric (Hermitian)
/// positive definite and `spd` says so, sparse LU otherwise (and as the fallback when Cholesky finds the matrix not positive definite).
/// The factorisation runs on the factorisation pool and cannot be interrupted (the stage line says so).
pub struct DirectSolver<T: FaerScalar> {
    fac: ShiftedFactorization<T>,
    n: usize,
}

impl<T: FaerScalar> DirectSolver<T> {
    pub fn new(a: &Csr<T>, spd: bool) -> Result<Self, String> {
        let n = a.nrows();
        if n == 0 || a.ncols() != n {
            return Err("the matrix to factorise must be square and not empty".into());
        }
        progress::check()?;
        let kind = if T::IS_COMPLEX { "complex " } else { "" };
        let _stage = progress::stage(format!("factorising a sparse {kind}{} matrix ({n} unknowns, {} non-zeros); this cannot be interrupted", if spd { "symmetric" } else { "general" }, a.nnz()));
        let mut fac = ShiftedFactorization::new(a, None, if spd { Definiteness::PositiveDefinite } else { Definiteness::General }).map_err(|e| format!("sparse factorisation failed: {e}"))?;
        fac.ordering = fill_ordering(n);
        crate::par::on_factor_pool(|| fac.set_shift(T::ZERO)).map_err(|e| format!("sparse factorisation failed (singular matrix?): {e}"))?;
        Ok(DirectSolver { fac, n })
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn solve(&self, b: &[T]) -> Vec<T> {
        assert_eq!(b.len(), self.n);
        let mut x = b.to_vec();
        self.fac.solve_block_in_place(&mut x, 1);
        x
    }

    /// Several right-hand sides; the independent triangular solves are spread over the threads of the shared pool.
    pub fn solve_many(&self, bs: &[Vec<T>]) -> Vec<Vec<T>> {
        let nt = if bs.len() > 1 && self.n > 5_000 { crate::par::threads().min(bs.len()) } else { 1 };
        if nt <= 1 {
            return bs.iter().map(|b| self.solve(b)).collect();
        }
        let chunk = bs.len().div_ceil(nt);
        let parts: Vec<Vec<Vec<T>>> = bs.par_chunks(chunk).map(|cols| cols.iter().map(|b| self.solve(b)).collect::<Vec<Vec<T>>>()).collect();
        parts.into_iter().flatten().collect()
    }
}

// ------------------------------------------------------------------------------------------------------------ eigenpairs

/// Options of the shift-and-invert eigensolvers. `tol` is the relative residual of the transformed problem
/// (`|OP x - nu x| / |nu x|`, SLEPc's test) that every returned pair meets.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EigenOptions {
    pub tol: f64,
    /// vectors per Krylov block; 0 picks `min(nev, 8)` (several vectors find every copy of a repeated eigenvalue)
    pub block: usize,
    /// largest Krylov basis; 0 picks the library default (about `2 nev`, at least `nev + 15` and room for the block)
    pub ncv: usize,
    /// solve the columns of a block on parallel threads (each thread runs the sequential triangular solve on a group of columns, as the
    /// old complex eigensolver did) when the problem has at least this many unknowns; 0 solves every block on one thread. The library's
    /// own parallel multi-RHS solve is not used (pathologically slow in faer 0.24.4, see `ThreadPolicy`).
    #[serde(default)]
    pub par_solve_min_n: usize,
}

impl Default for EigenOptions {
    fn default() -> Self {
        EigenOptions { tol: 1e-10, block: 0, ncv: 0, par_solve_min_n: 0 }
    }
}

/// Eigenpairs of `A x = lambda B x` nearest the shift, nearest first.
#[derive(Clone, Debug)]
pub struct Pairs<T: Scalar> {
    pub values: Vec<C64>,
    /// unit norm in the `B` inner product for a Hermitian pencil, unit 2-norm otherwise
    pub vectors: Vec<Vec<T>>,
    /// `|A x - lambda B x| / (|lambda| |x|)`, recomputed from the vectors
    pub residuals: Vec<f64>,
    /// the eigensolver's own relative error estimate of the transformed problem (`|OP x - nu x| / |nu x|`, what `EigenOptions::tol` bounds); the old stack reported this one
    pub errors: Vec<f64>,
    pub converged: bool,
    /// Krylov dimension used
    pub dim: usize,
    pub restarts: usize,
    pub matvecs: usize,
}

/// A projection applied to every operator result (removes a known null space of `A`). It must be a `B`-orthogonal projection.
pub type Projector<'a, T> = &'a dyn Fn(&mut Vec<T>);

/// Panic payload that carries a cancel out of the library's loops (the operator trait has no way to return an error).
struct CancelUnwind;

/// The shift-and-invert operator with a cancel check and an optional projection around it.
struct Guarded<'a, T: Scalar, O: LinearOperator<T>> {
    inner: O,
    project: Option<Projector<'a, T>>,
    _t: std::marker::PhantomData<T>,
}

impl<T: Scalar, O: LinearOperator<T>> Guarded<'_, T, O> {
    fn finish(&self, y: &mut [T], nrhs: usize) {
        let Some(p) = self.project else { return };
        let n = self.inner.dim();
        for c in 0..nrhs {
            let mut col = y[c * n..(c + 1) * n].to_vec();
            p(&mut col);
            y[c * n..(c + 1) * n].copy_from_slice(&col);
        }
    }
}

impl<T: Scalar, O: LinearOperator<T>> LinearOperator<T> for Guarded<'_, T, O> {
    fn dim(&self) -> usize {
        self.inner.dim()
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.apply_block(x, y, 1)
    }
    fn apply_block(&self, x: &[T], y: &mut [T], nrhs: usize) {
        if progress::check().is_err() {
            resume_unwind(Box::new(CancelUnwind));
        }
        self.inner.apply_block(x, y, nrhs);
        self.finish(y, nrhs);
    }
    fn apply_block_with_bx(&self, x: &[T], bx: &[T], y: &mut [T], nrhs: usize) {
        if progress::check().is_err() {
            resume_unwind(Box::new(CancelUnwind));
        }
        self.inner.apply_block_with_bx(x, bx, y, nrhs);
        self.finish(y, nrhs);
    }
}

/// The `nev` eigenpairs of `A x = lambda B x` nearest `sigma` (block Krylov-Schur on the shift-and-invert operator). `sym` says whether
/// the pencil is Hermitian with `B` positive definite and the shift real (then eigenvectors are `B`-orthonormal and the Cholesky
/// factorisation is available through `definiteness`), or general (complex symmetric PML pencils, tilted crystals). Cancels with the
/// error `"cancelled"`.
#[allow(clippy::too_many_arguments)]
pub fn eigs_near<T: FaerScalar>(a: &Csr<T>, b: &Csr<T>, sigma: T, nev: usize, sym: Symmetry, definiteness: Definiteness, opts: &EigenOptions, project: Option<Projector<'_, T>>) -> Result<Pairs<T>, String> {
    let n = a.nrows();
    if nev == 0 || nev >= n {
        return Err("number of eigenvalues must be between 1 and n-1".into());
    }
    if sym == Symmetry::General && !T::IS_COMPLEX {
        return Err("a general pencil needs complex arithmetic (use eigs_real_general_near)".into());
    }
    progress::check()?;
    let _stage = progress::stage(format!("shift-and-invert Krylov-Schur eigensolver ({nev} modes, {n} unknowns)"));
    let mut fac = ShiftedFactorization::new(a, Some(b), definiteness).map_err(|e| format!("sparse factorisation failed: {e}"))?;
    {
        let kind = if T::IS_COMPLEX { "complex " } else { "" };
        let _f = progress::stage(format!("factorising the {kind}shifted matrix ({n} unknowns); this cannot be interrupted"));
        crate::par::on_factor_pool(|| fac.set_shift(sigma)).map_err(|e| format!("sparse factorisation failed (singular matrix at the search shift?): {e}"))?;
    }
    progress::check()?;
    let block = if opts.block > 0 { opts.block } else { nev.min(8) };
    let kind = match sym {
        Symmetry::Hermitian => Kind::Hermitian,
        Symmetry::General => Kind::General,
    };
    let eps = EpsOptions {
        which: Which::LargestMagnitude,
        tol: opts.tol,
        block_size: block,
        ncv: if opts.ncv > 0 { Some(opts.ncv) } else { None },
        nev,
        ..match kind {
            Kind::Hermitian => EpsOptions::hermitian(nev),
            Kind::General => EpsOptions::general(nev),
        }
    };
    let op = Guarded { inner: ParShiftInvert { fac: &fac, b, min_n: opts.par_solve_min_n }, project, _t: std::marker::PhantomData };
    let weight: Option<&dyn LinearOperator<T>> = if sym == Symmetry::Hermitian { Some(b) } else { None };
    let res = match catch_unwind(AssertUnwindSafe(|| krylov_schur(&op, weight, None, &eps))) {
        Ok(r) => r.map_err(|e| format!("the eigensolver failed: {e}"))?,
        Err(payload) if payload.is::<CancelUnwind>() => return Err(progress::CANCELLED.into()),
        Err(payload) => resume_unwind(payload),
    };
    let s = sigma.to_c64();
    let values: Vec<C64> = res.values.iter().map(|&th| s + C64::ONE / th).collect();
    let (mut ax, mut bx) = (vec![T::ZERO; n], vec![T::ZERO; n]);
    let mut residuals = Vec::with_capacity(values.len());
    let mut vectors = Vec::with_capacity(values.len());
    for (i, &lam) in values.iter().enumerate() {
        let x = res.vectors.col(i);
        a.spmv(x, &mut ax);
        b.spmv(x, &mut bx);
        let r2: f64 = ax.iter().zip(&bx).map(|(p, q)| (p.to_c64() - lam * q.to_c64()).norm_sqr()).sum();
        residuals.push(r2.sqrt() / (lam.norm() * norm2(x)).max(1e-300));
        vectors.push(x.to_vec());
    }
    Ok(Pairs { values, vectors, residuals, errors: res.errors.iter().take(res.values.len()).copied().collect(), converged: res.reason == EpsReason::Converged, dim: res.ncv, restarts: res.iterations, matvecs: res.matvecs })
}

/// Real eigenpairs of a pencil with real, non-symmetric matrices (the 2D waveguide modes): the real eigenvalues nearest the real
/// shift, with real vectors. The library takes the pencil as complex; complex eigenvalues (conjugate pairs, not physical modes of the
/// formulation) are skipped, and the Krylov-Schur run asks for a few more pairs so that `nev` real ones remain.
pub fn eigs_real_general_near(a: &Csr<f64>, b: &Csr<f64>, sigma: f64, nev: usize, opts: &EigenOptions) -> Result<Pairs<f64>, String> {
    let n = a.nrows();
    if nev == 0 || nev >= n {
        return Err("number of eigenvalues must be between 1 and n-1".into());
    }
    let (ac, bc) = (a.to_complex(), b.to_complex());
    let mut ask = nev + 4;
    loop {
        let ask_now = ask.min(n - 1);
        let p = eigs_near(&ac, &bc, C64::new(sigma, 0.0), ask_now, Symmetry::General, Definiteness::General, opts, None)?;
        if std::env::var_os("NYX_DEBUG_EIG").is_some() { eprintln!("sigma {sigma:e} ask {ask_now} conv {} values {:?}", p.converged, p.values.iter().take(8).collect::<Vec<_>>()); }
        let mut out = Pairs { values: vec![], vectors: vec![], residuals: vec![], errors: vec![], converged: p.converged, dim: p.dim, restarts: p.restarts, matvecs: p.matvecs };
        for (i, &lam) in p.values.iter().enumerate() {
            if lam.im.abs() > 1e-8 * lam.re.abs().max(1e-300) {
                continue;
            }
            // the real eigenvector: rotate the phase so that the largest entry is real, then drop the (tiny) imaginary part
            let big = p.vectors[i].iter().copied().fold(C64::new(0.0, 0.0), |m, z| if z.norm_sqr() > m.norm_sqr() { z } else { m });
            let rot = big.conj() / big.norm();
            let mut x: Vec<f64> = p.vectors[i].iter().map(|z| (z * rot).re).collect();
            let xn = norm2(&x);
            x.iter_mut().for_each(|t| *t /= xn);
            out.values.push(C64::new(lam.re, 0.0));
            out.vectors.push(x);
            out.residuals.push(p.residuals[i]);
            out.errors.push(p.errors.get(i).copied().unwrap_or(0.0));
            if out.values.len() == nev {
                break;
            }
        }
        if out.values.len() == nev {
            return Ok(out);
        }
        if ask_now >= n - 1 || ask > 8 * nev + 32 {
            out.converged = false;
            return if out.values.is_empty() { Err("no real eigenvalue was found near the shift".into()) } else { Ok(out) };
        }
        ask = 2 * ask;
    }
}

// ------------------------------------------------------------------------------------------------------------ small dense problems

/// Dense matrix as row-major `Vec<f64>` with `n` columns, as a faer matrix.
fn to_faer(a: &[f64], nrows: usize, ncols: usize) -> faer::Mat<f64> {
    faer::Mat::from_fn(nrows, ncols, |i, j| a[i * ncols + j])
}

/// Eigen decomposition of a symmetric dense matrix (row-major, only the lower triangle is read): ascending eigenvalues and the
/// eigenvectors (`vecs[j]` is the j-th).
pub fn dense_eigh(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let m = to_faer(a, n, n);
    let e = m.self_adjoint_eigen(faer::Side::Lower).map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    let (u, s) = (e.U(), e.S());
    let vals: Vec<f64> = (0..n).map(|i| s.column_vector()[i]).collect();
    let vecs: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| u[(i, j)]).collect()).collect();
    Ok((vals, vecs))
}

/// Generalised symmetric-definite problem `K x = lambda M x` (dense, row-major, `M` positive definite) through the Cholesky factor
/// of `M`: ascending eigenvalues and `M`-orthonormal eigenvectors.
pub fn dense_gen_eigh(k: &[f64], m: &[f64], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let mm = to_faer(m, n, n);
    let llt = mm.llt(faer::Side::Lower).map_err(|e| format!("mass matrix is not positive definite: {e:?}"))?;
    let l = llt.L();
    // C = L^-1 K L^-T
    let mut c = to_faer(k, n, n);
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, c.as_mut(), faer::get_global_parallelism());
    let mut ct = c.transpose().to_owned();
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, ct.as_mut(), faer::get_global_parallelism());
    // symmetrise against round-off
    let cs = faer::Mat::from_fn(n, n, |i, j| 0.5 * (ct[(i, j)] + ct[(j, i)]));
    let e = cs.self_adjoint_eigen(faer::Side::Lower).map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    let (u, s) = (e.U(), e.S());
    let vals: Vec<f64> = (0..n).map(|i| s.column_vector()[i]).collect();
    // x = L^-T y
    let mut y = faer::Mat::from_fn(n, n, |i, j| u[(i, j)]);
    faer::linalg::triangular_solve::solve_upper_triangular_in_place(l.transpose(), y.as_mut(), faer::get_global_parallelism());
    let vecs: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| y[(i, j)]).collect()).collect();
    Ok((vals, vecs))
}

// ------------------------------------------------------------------------------------------------------------ cancel helpers

/// Call inside an operator or preconditioner application: leaves the library's loop with a cancel when the job was cancelled.
pub fn cancel_point() {
    if progress::check().is_err() {
        resume_unwind(Box::new(CancelUnwind));
    }
}

/// Runs a library solve, turning a cancel raised by [`cancel_point`] into `Err("cancelled")`.
pub fn catching_cancel<R>(f: impl FnOnce() -> R) -> Result<R, String> {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(r) => Ok(r),
        Err(payload) if payload.is::<CancelUnwind>() => Err(progress::CANCELLED.into()),
        Err(payload) => resume_unwind(payload),
    }
}

// ------------------------------------------------------------------------------------------------------------ kernel-facing wrappers

pub use crate::solvers::arnoldi_real::{EigOpts, EigPairs};
use crate::solvers::csr::{CCsr, RCsr};

/// Start of the error text returned when the LOBPCG stall guard gives up (the caller falls back to a factorisation).
pub const STALLED: &str = "LOBPCG stalled";

/// Eigensolver tolerance override for experiments on the tolerance rule (gate G2): `NYX_EIG_TOL=1e-12` replaces the tolerance every
/// kernel asks for. Not set: each kernel keeps its own (1e-9 to 1e-11).
fn tol_override() -> Option<f64> {
    use std::sync::OnceLock;
    static T: OnceLock<Option<f64>> = OnceLock::new();
    *T.get_or_init(|| std::env::var("NYX_EIG_TOL").ok().and_then(|v| v.trim().parse::<f64>().ok()).filter(|t| *t > 0.0))
}

fn opts_of(tol: f64, block: usize) -> EigenOptions {
    EigenOptions { tol: tol_override().unwrap_or(tol), block, ncv: 0, par_solve_min_n: 0 }
}

/// Column-parallel solves for the complex stages: the old complex eigensolver spread the columns of a block over the pool for n > 5000. `NYX_PAR_SOLVE_MIN_N=<n>` changes
/// the size from which it happens, `0` turns it off (sequential solves), for benchmarking.
fn par_solve_min_n_complex() -> usize {
    use std::sync::OnceLock;
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| std::env::var("NYX_PAR_SOLVE_MIN_N").ok().and_then(|v| v.trim().parse::<usize>().ok()).unwrap_or(5000))
}

fn definiteness(spd: bool) -> Definiteness {
    if spd { Definiteness::PositiveDefinite } else { Definiteness::General }
}

fn real_pairs(p: Pairs<f64>) -> EigPairs {
    EigPairs { values: p.values.iter().map(|z| z.re).collect(), vectors: p.vectors, dim: p.dim, converged: p.converged, complex_found: false }
}

/// Real symmetric pencil `A x = lambda M x` with `M` positive definite, eigenvalues nearest `sigma`; a shift at or below zero is taken to
/// leave `A - sigma M` positive definite (Cholesky), otherwise LU is used. `project` removes a known null space from every operator result.
pub fn eigs_sym(a: &RCsr, m: &RCsr, sigma: f64, nev: usize, opts: &EigOpts, project: Option<Projector<f64>>) -> Result<EigPairs, String> {
    let p = eigs_near(a, m, sigma, nev, Symmetry::Hermitian, definiteness(sigma <= 0.0), &opts_of(opts.tol, opts.block), project)?;
    Ok(real_pairs(p))
}

/// Like [`eigs_sym`] for an indefinite `A - sigma M` (always LU).
pub fn eigs_sym_indefinite(a: &RCsr, m: &RCsr, sigma: f64, nev: usize, opts: &EigOpts, project: Option<Projector<f64>>) -> Result<EigPairs, String> {
    let p = eigs_near(a, m, sigma, nev, Symmetry::Hermitian, Definiteness::General, &opts_of(opts.tol, opts.block), project)?;
    Ok(real_pairs(p))
}

/// Real eigenvalues and real vectors of a real non-symmetric pencil (2D waveguide modes), in real arithmetic
/// (see `arnoldi_real` for why the complex Krylov-Schur is not used here).
pub fn eigs_general(a: &RCsr, b: &RCsr, sigma: f64, nev: usize, opts: &EigOpts) -> Result<EigPairs, String> {
    crate::solvers::arnoldi_real::eigs_general(a, b, sigma, nev, opts)
}

/// Eigenpairs of a complex pencil.
#[derive(Clone, Debug)]
pub struct CEigPairs {
    pub values: Vec<C64>,
    pub vectors: Vec<Vec<C64>>,
    /// the eigensolver's relative error estimate of the transformed problem (what the tolerance bounds)
    pub residuals: Vec<f64>,
    pub dim: usize,
    pub converged: bool,
}

/// Complex pencil (PML, tilted crystals, Floquet phases), nearest the complex shift; `project` as in [`eigs_sym`].
pub fn eigs_complex_projected(a: &CCsr, b: &CCsr, sigma: C64, nev: usize, tol: f64, _max_dim: usize, project: Option<Projector<C64>>) -> Result<CEigPairs, String> {
    let opts = EigenOptions { par_solve_min_n: par_solve_min_n_complex(), ..opts_of(tol, 0) };
    let p = eigs_near(a, b, sigma, nev, Symmetry::General, Definiteness::General, &opts, project)?;
    Ok(CEigPairs { values: p.values, vectors: p.vectors, residuals: p.errors, dim: p.dim, converged: p.converged })
}

pub fn eigs_complex(a: &CCsr, b: &CCsr, sigma: C64, nev: usize, tol: f64, max_dim: usize) -> Result<CEigPairs, String> {
    eigs_complex_projected(a, b, sigma, nev, tol, max_dim, None)
}

/// A preconditioner given as a closure over free-dof vectors (one multigrid cycle).
struct ClosurePc<'a>(&'a dyn Fn(&[f64]) -> Vec<f64>);

impl crate::solvers::linop::Preconditioner<f64> for ClosurePc<'_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        cancel_point();
        z.copy_from_slice(&(self.0)(r));
    }
}

/// The `nev` lowest eigenpairs of `K x = lambda M x` by LOBPCG with the given preconditioner. With `stall_guard` the run gives up (error
/// starting with [`STALLED`]) as soon as the worst residual improves by less than 10 % over 10 iterations, so that the caller can fall
/// back to a factorisation.
#[allow(clippy::too_many_arguments)]
pub fn lobpcg_guarded(k: &RCsr, m: &RCsr, nev: usize, extra: usize, precond: &dyn Fn(&[f64]) -> Vec<f64>, tol: f64, max_iter: usize, stall_guard: bool) -> Result<EigPairs, String> {
    use crate::solvers::lobpcg::{LobpcgOptions, LobpcgReason, lobpcg};
    let n = k.nrows();
    let _stage = progress::stage(format!("LOBPCG eigensolver with multigrid preconditioner ({nev} modes, {n} unknowns)"));
    let opts = LobpcgOptions { block: nev + extra, tol, max_it: max_iter, stall_window: if stall_guard { 10 } else { max_iter.max(1) }, stall_gain: if stall_guard { 0.1 } else { 0.01 }, ..LobpcgOptions::new(nev) };
    let pc = ClosurePc(precond);
    let res = catching_cancel(|| lobpcg(k, Some(m), Some(&pc), None, &[], &opts))?.map_err(|e| format!("the LOBPCG eigensolver failed: {e}"))?;
    if res.reason == LobpcgReason::Stalled {
        return Err(format!("{STALLED} after {} iterations", res.iterations));
    }
    Ok(EigPairs { values: res.values, vectors: res.vectors, dim: nev + extra, converged: res.reason == LobpcgReason::Converged, complex_found: false })
}

// ------------------------------------------------------------------------------------------------------------ column-parallel solves

/// `x -> (A - sigma B)^-1 B x` like [`ShiftInvertOperator`], with the columns of a block solved on parallel threads when
/// `n >= min_n` (each thread runs faer's sequential triangular solve on its group of columns: independent solves against one shared
/// factorisation). `min_n == 0` or a single column is the library operator unchanged.
struct ParShiftInvert<'a, T: FaerScalar> {
    fac: &'a ShiftedFactorization<T>,
    b: &'a Csr<T>,
    min_n: usize,
}

impl<T: FaerScalar> ParShiftInvert<'_, T> {
    fn solve_block(&self, y: &mut [T], nrhs: usize) {
        let n = self.fac.dim();
        let nt = if self.min_n > 0 && nrhs > 1 && n >= self.min_n { crate::par::threads().min(nrhs) } else { 1 };
        if nt <= 1 {
            self.fac.solve_block_in_place(&mut y[..n * nrhs], nrhs);
            return;
        }
        let cols = nrhs.div_ceil(nt);
        y[..n * nrhs].par_chunks_mut(cols * n).for_each(|chunk| {
            let k = chunk.len() / n;
            self.fac.solve_block_in_place(chunk, k);
        });
    }
}

impl<T: FaerScalar> LinearOperator<T> for ParShiftInvert<'_, T> {
    fn dim(&self) -> usize {
        self.fac.dim()
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.apply_block(x, y, 1)
    }
    fn apply_block(&self, x: &[T], y: &mut [T], nrhs: usize) {
        let n = self.fac.dim();
        for c in 0..nrhs {
            self.b.spmv(&x[c * n..(c + 1) * n], &mut y[c * n..(c + 1) * n]);
        }
        self.solve_block(y, nrhs);
    }
    fn apply_block_with_bx(&self, _x: &[T], bx: &[T], y: &mut [T], nrhs: usize) {
        let n = self.fac.dim();
        y[..n * nrhs].copy_from_slice(&bx[..n * nrhs]);
        self.solve_block(y, nrhs);
    }
}

/// Fill-reducing ordering of the Cholesky path: AMD, as the old stack. The library's `Auto` (also tries nested dissection, keeps the one with
/// fewer factor non-zeros) was measured in phase 3 and lost: the 205 479-unknown 3D direct solve took 13.5 and 17.8 s with `Auto` against 4.7
/// and 10.9 s with AMD (two runs each, loaded machine), and the small 2D solves 60 to 80 % longer (the second symbolic analysis). The
/// library's 30 % less fill for 3D Laplacians does not carry over to P2 elasticity with the parallel factorisation.
/// `NYX_ORDERING=amd|auto|nd` overrides the choice (experiments).
fn fill_ordering(_n: usize) -> crate::solvers::ordering::FillOrdering {
    use crate::solvers::ordering::FillOrdering;
    match std::env::var("NYX_ORDERING").ok().as_deref() {
        Some("amd") => FillOrdering::Amd,
        Some("auto") => FillOrdering::Auto,
        Some("nd") => FillOrdering::NestedDissection,
        _ => FillOrdering::Amd,
    }
}
