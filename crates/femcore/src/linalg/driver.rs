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
use crate::linalg::csr::Csr;
pub use crate::linalg::direct::{FactorBackend, FactorStats};
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::linalg::eps::{EpsOptions, EpsReason, Kind, Which, krylov_schur};
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::linop::LinearOperator;
use crate::linalg::problem::Symmetry;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::norm2;
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

    /// Replaces the numeric values of `a` while reusing the cached symbolic analysis
    /// (requiring an identical sparsity pattern), and re-runs numeric factorisation on the factorisation pool.
    pub fn update_matrix(&mut self, a: &Csr<T>) -> Result<(), String> {
        if a.nrows() != self.n || a.ncols() != self.n {
            return Err("updated matrix dimension does not match the solver dimension".into());
        }
        progress::check()?;
        self.fac
            .replace_matrices(a, None)
            .map_err(|e| format!("sparse matrix update failed: {e}"))?;
        crate::par::on_factor_pool(|| self.fac.set_shift(T::ZERO))
            .map_err(|e| format!("sparse factorisation failed (singular matrix?): {e}"))?;
        Ok(())
    }

    /// Factorisation statistics (including `symbolic_count` and `numeric_count`).
    pub fn stats(&self) -> FactorStats {
        self.fac.stats()
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
    fac.ordering = fill_ordering(n);
    {
        let kind = if T::IS_COMPLEX { "complex " } else { "" };
        let _f = progress::stage(format!("factorising the {kind}shifted matrix ({n} unknowns); this cannot be interrupted"));
        crate::par::on_factor_pool(|| fac.set_shift(sigma)).map_err(|e| format!("sparse factorisation failed (singular matrix at the search shift?): {e}"))?;
    }
    progress::check()?;
    let block = if opts.block > 0 { opts.block } else { nev.min(8) };
    // a warm start (P5b): the modes of the previous solve of a sweep are the start block, all of them, so the block holds at least `nev` columns
    let ticket = crate::warm::begin("shift_invert", n, nev);
    let warm_start: Option<Vec<T>> = ticket.as_ref().and_then(|t| t.start_columns::<T>());
    let block = if warm_start.is_some() { block.max(nev) } else { block };
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
    let res = match catch_unwind(AssertUnwindSafe(|| krylov_schur(&op, weight, warm_start.as_deref(), &eps))) {
        Ok(r) => r.map_err(|e| format!("the eigensolver failed: {e}"))?,
        Err(payload) if payload.is::<CancelUnwind>() => return Err(progress::CANCELLED.into()),
        Err(payload) => resume_unwind(payload),
    };
    if let Some(t) = ticket {
        if res.reason == EpsReason::Converged {
            t.finish((0..res.values.len()).map(|i| res.vectors.col(i)), res.iterations);
        }
    }
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
        ask *= 2;
    }
}

// ------------------------------------------------------------------------------------------------------------ small dense problems

/// Dense matrix as row-major `Vec<f64>` with `n` columns, as a faer matrix.
fn to_faer(a: &[f64], nrows: usize, ncols: usize) -> faer::Mat<f64> {
    faer::Mat::from_fn(nrows, ncols, |i, j| a[i * ncols + j])
}

/// Eigen decomposition of a symmetric dense matrix (row-major, only the lower triangle is read): ascending eigenvalues and the
/// eigenvectors (`vecs[j]` is the j-th). Uses `Par::Seq` for `n <= 512` (`ADR 11`).
pub fn dense_eigh(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::evd::{ComputeEigenvectors, self_adjoint_evd, self_adjoint_evd_scratch};

    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let m = faer::Mat::from_fn(n, n, |i, j| if i >= j { a[i * n + j] } else { a[j * n + i] });
    let par = crate::linalg::faer_glue::dense_par_dim(n);
    let mut u = faer::Mat::<f64>::zeros(n, n);
    let mut s = faer::diag::Diag::<f64>::zeros(n);
    let mut mem = MemBuffer::new(self_adjoint_evd_scratch::<f64>(
        n,
        ComputeEigenvectors::Yes,
        par,
        Default::default(),
    ));
    let res = self_adjoint_evd(
        m.as_ref(),
        s.as_mut(),
        Some(u.as_mut()),
        par,
        MemStack::new(&mut mem),
        Default::default(),
    );
    crate::linalg::faer_glue::clear_upper_vector_state();
    res.map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    let vals: Vec<f64> = (0..n).map(|i| s[i]).collect();
    let vecs: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| u[(i, j)]).collect()).collect();
    Ok((vals, vecs))
}

/// Generalised symmetric-definite problem `K x = lambda M x` (dense, row-major, `M` positive definite) through the Cholesky factor
/// of `M`: ascending eigenvalues and `M`-orthonormal eigenvectors. Uses `Par::Seq` for `n <= 512` (`ADR 11`).
pub fn dense_gen_eigh(k: &[f64], m: &[f64], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::cholesky::llt::factor::{cholesky_in_place, cholesky_in_place_scratch};
    use faer::linalg::evd::{ComputeEigenvectors, self_adjoint_evd, self_adjoint_evd_scratch};

    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let par = crate::linalg::faer_glue::dense_par_dim(n);
    let mut mm = to_faer(m, n, n);
    {
        let mut llt_mem = MemBuffer::new(cholesky_in_place_scratch::<f64>(n, par, Default::default()));
        let res = cholesky_in_place(
            mm.as_mut(),
            Default::default(),
            par,
            MemStack::new(&mut llt_mem),
            Default::default(),
        );
        crate::linalg::faer_glue::clear_upper_vector_state();
        res.map_err(|e| format!("mass matrix is not positive definite: {e:?}"))?;
    }
    let l = mm.as_ref();
    // C = L^-1 K L^-T
    let mut c = to_faer(k, n, n);
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, c.as_mut(), par);
    let mut ct = c.transpose().to_owned();
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, ct.as_mut(), par);
    // symmetrise against round-off
    let cs = faer::Mat::from_fn(n, n, |i, j| 0.5 * (ct[(i, j)] + ct[(j, i)]));
    let mut u = faer::Mat::<f64>::zeros(n, n);
    let mut s = faer::diag::Diag::<f64>::zeros(n);
    {
        let mut evd_mem = MemBuffer::new(self_adjoint_evd_scratch::<f64>(
            n,
            ComputeEigenvectors::Yes,
            par,
            Default::default(),
        ));
        let res = self_adjoint_evd(
            cs.as_ref(),
            s.as_mut(),
            Some(u.as_mut()),
            par,
            MemStack::new(&mut evd_mem),
            Default::default(),
        );
        crate::linalg::faer_glue::clear_upper_vector_state();
        res.map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    }
    let vals: Vec<f64> = (0..n).map(|i| s[i]).collect();
    // x = L^-T y
    faer::linalg::triangular_solve::solve_upper_triangular_in_place(l.transpose(), u.as_mut(), par);
    crate::linalg::faer_glue::clear_upper_vector_state();
    let vecs: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| u[(i, j)]).collect()).collect();
    Ok((vals, vecs))
}

/// Thin SVD of a complex matrix given by its columns (`cols[j]` has length `n`):
/// returns `(singular_values, left_singular_vectors)` in decreasing singular-value order.
pub fn thin_svd_c64(cols: &[Vec<C64>]) -> Result<(Vec<f64>, Vec<Vec<C64>>), String> {
    if cols.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let n = cols[0].len();
    if cols.iter().any(|c| c.len() != n) {
        return Err("thin_svd_c64: all columns must have the same length".into());
    }
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let k = cols.len();
    if n >= k {
        let mut data = Vec::with_capacity(n * k);
        for col in cols {
            data.extend_from_slice(col);
        }
        let (s, u_flat) =
            crate::linalg::faer_glue::thin_svd_c64(&data, n, k).map_err(|e| format!("thin SVD failed: {e}"))?;
        let u_cols = (0..k).map(|j| u_flat[j * n..(j + 1) * n].to_vec()).collect();
        Ok((s, u_cols))
    } else {
        if cols.iter().flatten().any(|v| !v.is_finite()) {
            return Err("thin SVD failed: matrix contains NaN or infinity".into());
        }
        use faer::dyn_stack::{MemBuffer, MemStack};
        use faer::linalg::svd::{ComputeSvdVectors, svd, svd_scratch};
        let m = faer::Mat::<faer::c64>::from_fn(n, k, |i, j| cols[j][i]);
        let par = crate::linalg::faer_glue::dense_par(n * k);
        let mut s_diag = faer::diag::Diag::<faer::c64>::zeros(n);
        let mut u = faer::Mat::<faer::c64>::zeros(n, n);
        let mut mem = MemBuffer::new(svd_scratch::<faer::c64>(
            n,
            k,
            ComputeSvdVectors::Thin,
            ComputeSvdVectors::No,
            par,
            Default::default(),
        ));
        let res = svd(
            m.as_ref(),
            s_diag.as_mut(),
            Some(u.as_mut()),
            None,
            par,
            MemStack::new(&mut mem),
            Default::default(),
        )
        .map_err(|_| "thin SVD failed to converge".to_string());
        crate::linalg::faer_glue::clear_upper_vector_state();
        res?;
        let s: Vec<f64> = (0..n).map(|i| s_diag[i].re).collect();
        let u_cols: Vec<Vec<C64>> = (0..n).map(|j| (0..n).map(|i| u[(i, j)]).collect()).collect();
        Ok((s, u_cols))
    }
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

pub use crate::linalg::arnoldi_real::{EigOpts, EigPairs};
use crate::linalg::csr::{CCsr, RCsr};

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
    crate::linalg::arnoldi_real::eigs_general(a, b, sigma, nev, opts)
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

/// Complex quadratic eigenproblem `(A0 + lambda A1 + lambda^2 A2) u = 0`, returning up to `nev` eigenpairs nearest `sigma`.
///
/// Linearises the shifted pencil `P(sigma + mu) u = (P(sigma) + mu P'(sigma) + mu^2 A2) u = 0` (`lambda = sigma + 1/theta`)
/// through the `2n x 2n` shift-and-invert companion operator `T [u; w] = [w; -P(sigma)^-1 (A2 u + P'(sigma) w)]`, factorising
/// only the `n x n` matrix `P(sigma) = A0 + sigma A1 + sigma^2 A2` once (`A2` may be singular; zero-`theta` roots at infinity
/// sit at the origin of `T` and are skipped).
pub fn eigs_quadratic_complex(a0: &CCsr, a1: &CCsr, a2: &CCsr, sigma: C64, nev: usize, tol: f64) -> Result<CEigPairs, String> {
    let n = a0.nrows();
    if n == 0 || a0.ncols() != n || a1.nrows() != n || a1.ncols() != n || a2.nrows() != n || a2.ncols() != n {
        return Err("quadratic eigenproblem matrices must be square and of the same non-zero dimension".into());
    }
    if nev == 0 || nev >= 2 * n {
        return Err("number of eigenvalues must be between 1 and 2n-1".into());
    }
    progress::check()?;
    let _stage = progress::stage(format!("quadratic shift-and-invert Krylov-Schur eigensolver ({nev} modes, {n} unknowns)"));
    let p_sigma = a0.add_scaled(a1, sigma).add_scaled(a2, sigma * sigma);
    let p_prime = a1.add_scaled(a2, C64::new(2.0, 0.0) * sigma);
    let mut fac = ShiftedFactorization::new(&p_sigma, None, Definiteness::General)
        .map_err(|e| format!("sparse factorisation failed: {e}"))?;
    fac.ordering = fill_ordering(n);
    {
        let _f = progress::stage(format!(
            "factorising the complex quadratic shifted matrix ({n} unknowns); this cannot be interrupted"
        ));
        crate::par::on_factor_pool(|| fac.set_shift(C64::ZERO))
            .map_err(|e| format!("sparse factorisation failed (singular quadratic pencil at the search shift?): {e}"))?;
    }
    progress::check()?;

    let tol_eff = tol_override().unwrap_or(tol);
    let mut block = nev.min(8);
    let ticket = crate::warm::begin("qep_shift_invert", 2 * n, nev);
    let warm_start: Option<Vec<C64>> = ticket.as_ref().and_then(|t| t.start_columns::<C64>());
    if warm_start.is_some() {
        block = block.max(nev);
    }
    let eps = EpsOptions {
        which: Which::LargestMagnitude,
        tol: tol_eff,
        block_size: block,
        ncv: None,
        nev,
        ..EpsOptions::general(nev)
    };
    let op = Guarded {
        inner: QepShiftInvert { fac: &fac, a2, p_prime: &p_prime, min_n: par_solve_min_n_complex() },
        project: None,
        _t: std::marker::PhantomData,
    };
    let res = match catch_unwind(AssertUnwindSafe(|| krylov_schur(&op, None, warm_start.as_deref(), &eps))) {
        Ok(r) => r.map_err(|e| format!("the quadratic eigensolver failed: {e}"))?,
        Err(payload) if payload.is::<CancelUnwind>() => return Err(progress::CANCELLED.into()),
        Err(payload) => resume_unwind(payload),
    };
    if let Some(t) = ticket {
        if res.reason == EpsReason::Converged {
            t.finish((0..res.values.len()).map(|i| res.vectors.col(i)), res.iterations);
        }
    }

    let (mut r0, mut r1, mut r2) = (vec![C64::ZERO; n], vec![C64::ZERO; n], vec![C64::ZERO; n]);
    let mut pairs: Vec<(C64, Vec<C64>, f64)> = Vec::with_capacity(res.values.len());
    for (i, &theta) in res.values.iter().enumerate() {
        if theta.norm() < 1e-15 {
            continue;
        }
        let lam = sigma + C64::ONE / theta;
        let col = res.vectors.col(i);
        let u_raw = &col[..n];
        let unorm = norm2(u_raw);
        if unorm <= 1e-300 {
            continue;
        }
        let u: Vec<C64> = u_raw.iter().map(|&z| z / unorm).collect();
        a0.spmv(&u, &mut r0);
        a1.spmv(&u, &mut r1);
        a2.spmv(&u, &mut r2);
        let lam2 = lam * lam;
        let num = (0..n).map(|k| (r0[k] + lam * r1[k] + lam2 * r2[k]).norm_sqr()).sum::<f64>().sqrt();
        let den = norm2(&r0) + lam.norm() * norm2(&r1) + lam2.norm() * norm2(&r2) + 1e-300;
        pairs.push((lam, u, num / den));
    }
    pairs.sort_by(|a, b| (a.0 - sigma).norm().total_cmp(&(b.0 - sigma).norm()));
    let mut values = Vec::with_capacity(pairs.len());
    let mut vectors = Vec::with_capacity(pairs.len());
    let mut residuals = Vec::with_capacity(pairs.len());
    for (v, u, r) in pairs {
        values.push(v);
        vectors.push(u);
        residuals.push(r);
    }
    Ok(CEigPairs { values, vectors, residuals, dim: res.ncv, converged: res.reason == EpsReason::Converged })
}

/// A preconditioner given as a closure over free-dof vectors (one multigrid cycle).
struct ClosurePc<'a>(&'a dyn Fn(&[f64]) -> Vec<f64>);

impl crate::linalg::linop::Preconditioner<f64> for ClosurePc<'_> {
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
    use crate::linalg::lobpcg::{LobpcgOptions, LobpcgReason, lobpcg};
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
        if self.min_n == 0 || nrhs <= 1 || n < self.min_n {
            self.fac.solve_block_in_place(&mut y[..n * nrhs], nrhs);
            return;
        }
        // Always solve 1 column per call when column-parallel mode is enabled so `NYX_THREADS=1`
        // and `NYX_THREADS=N` execute the exact same triangular-solve kernel (`ADR 11`).
        if crate::par::threads() <= 1 {
            for chunk in y[..n * nrhs].chunks_mut(n) {
                self.fac.solve_block_in_place(chunk, 1);
            }
        } else {
            y[..n * nrhs].par_chunks_mut(n).for_each(|chunk| {
                self.fac.solve_block_in_place(chunk, 1);
            });
        }
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

/// Companion shift-and-invert operator `T [u; w] = [w; -P(sigma)^-1 (A2 u + P'(sigma) w)]` on `C^{2n}` for
/// [`eigs_quadratic_complex`].
struct QepShiftInvert<'a> {
    fac: &'a ShiftedFactorization<C64>,
    a2: &'a CCsr,
    p_prime: &'a CCsr,
    min_n: usize,
}

impl QepShiftInvert<'_> {
    fn solve_block(&self, rhs: &mut [C64], nrhs: usize) {
        let n = self.fac.dim();
        if self.min_n == 0 || nrhs <= 1 || n < self.min_n {
            self.fac.solve_block_in_place(&mut rhs[..n * nrhs], nrhs);
            return;
        }
        // Always solve 1 column per call when column-parallel mode is enabled so `NYX_THREADS=1`
        // and `NYX_THREADS=N` execute the exact same triangular-solve kernel (`ADR 11`).
        if crate::par::threads() <= 1 {
            for chunk in rhs[..n * nrhs].chunks_mut(n) {
                self.fac.solve_block_in_place(chunk, 1);
            }
        } else {
            rhs[..n * nrhs].par_chunks_mut(n).for_each(|chunk| {
                self.fac.solve_block_in_place(chunk, 1);
            });
        }
    }
}

impl LinearOperator<C64> for QepShiftInvert<'_> {
    fn dim(&self) -> usize {
        2 * self.fac.dim()
    }
    fn apply(&self, x: &[C64], y: &mut [C64]) {
        self.apply_block(x, y, 1)
    }
    fn apply_block(&self, x: &[C64], y: &mut [C64], nrhs: usize) {
        let n = self.fac.dim();
        let n2 = 2 * n;
        let mut rhs = vec![C64::ZERO; n * nrhs];
        let mut tmp = vec![C64::ZERO; n];
        for c in 0..nrhs {
            let xc = &x[c * n2..(c + 1) * n2];
            let (u, w) = (&xc[..n], &xc[n..n2]);
            y[c * n2..c * n2 + n].copy_from_slice(w);
            let rc = &mut rhs[c * n..(c + 1) * n];
            self.a2.spmv(u, rc);
            self.p_prime.spmv(w, &mut tmp);
            for i in 0..n {
                rc[i] = -(rc[i] + tmp[i]);
            }
        }
        self.solve_block(&mut rhs, nrhs);
        for c in 0..nrhs {
            y[c * n2 + n..(c + 1) * n2].copy_from_slice(&rhs[c * n..(c + 1) * n]);
        }
    }
}

// ------------------------------------------------------------------------------------------------------------ P2 Nédélec iterative curl-curl solvers (FR-16 Pt 2)

pub use crate::linalg::ams::{AmsP2, AmsP2Options, AmsP2Topology};

/// Options for [`solve_curl_curl_p2_iterative`] (FGMRES/GMRES preconditioned by [`AmsP2`]).
#[derive(Clone, Debug)]
pub struct IterativeCurlCurlOptions {
    /// Target relative residual `‖A x − b‖₂ / ‖b‖₂`.
    pub rtol: f64,
    /// Absolute residual stopping tolerance.
    pub atol: f64,
    /// Maximum Krylov iterations.
    pub max_it: usize,
    /// Restart dimension `m` of FGMRES/GMRES(m).
    pub restart: usize,
    /// Use Flexible GMRES (`true`, default) instead of standard right-preconditioned GMRES.
    pub flexible: bool,
    /// Options of the $P_2$ Hiptmair–Xu auxiliary-space preconditioner.
    pub ams: AmsP2Options,
}

impl Default for IterativeCurlCurlOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-9,
            atol: 1e-14,
            max_it: 400,
            restart: 40,
            flexible: true,
            ams: AmsP2Options::default(),
        }
    }
}

/// Output of [`solve_curl_curl_p2_iterative`].
#[derive(Clone, Debug)]
pub struct IterativeSolveReport<T: Scalar> {
    pub x: Vec<T>,
    pub iterations: usize,
    /// Independently recomputed true relative residual `‖A x − b‖₂ / ‖b‖₂`.
    pub relative_residual: f64,
    pub converged: bool,
    /// Total sparse non-zeros stored by the `AmsP2` preconditioner hierarchy.
    pub preconditioner_nnz: usize,
}

struct CancelPc<'a, T: Scalar, P: crate::linalg::linop::Preconditioner<T>>(&'a P, std::marker::PhantomData<T>);

impl<T: Scalar, P: crate::linalg::linop::Preconditioner<T>> crate::linalg::linop::Preconditioner<T>
    for CancelPc<'_, T, P>
{
    fn apply(&self, r: &[T], z: &mut [T]) {
        cancel_point();
        self.0.apply(r, z);
    }
}

/// Solves a 3D $P_2$ Nédélec curl-curl linear system `A x = b` (driven frequency-domain or
/// shifted Helmholtz/Maxwell) iteratively with FGMRES/GMRES preconditioned by [`AmsP2`].
///
/// When `A = K − σ M` is indefinite, pass a positive-definite or complex-shifted auxiliary
/// matrix as `pc_matrix` (e.g. `Some(&(K + |σ| M))` or `Some(&(K − (σ + iγ) M))`); when `A` is
/// already SPD or complex-damped, pass `pc_matrix = None` to build `AmsP2` directly on `A`.
pub fn solve_curl_curl_p2_iterative<T: FaerScalar>(
    a: &Csr<T>,
    pc_matrix: Option<&Csr<T>>,
    topo: &AmsP2Topology,
    rhs: &[T],
    opts: &IterativeCurlCurlOptions,
) -> Result<IterativeSolveReport<T>, String> {
    use crate::linalg::ksp::{Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances};
    let n = a.nrows();
    if n == 0 || a.ncols() != n || rhs.len() != n {
        return Err("solve_curl_curl_p2_iterative: matrix must be square and match rhs length".into());
    }
    progress::check()?;
    let _stage = progress::stage(format!("iterative P2 curl-curl FGMRES+AMS ({n} unknowns)"));
    let pc_mat = pc_matrix.unwrap_or(a);
    let ams = AmsP2::new(pc_mat, topo, &opts.ams).map_err(|e| format!("AmsP2 setup failed: {e}"))?;
    let preconditioner_nnz = ams.nnz_storage();
    let gmres_opts = GmresOptions {
        tol: Tolerances { rtol: opts.rtol, atol: opts.atol, max_it: opts.max_it, ..Default::default() },
        restart: opts.restart.max(2),
        side: PcSide::Right,
        orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
        flexible: opts.flexible,
        ..Default::default()
    };
    let mut solver = Gmres::<T>::new(gmres_opts);
    let mut x = vec![T::ZERO; n];
    let pc = CancelPc(&ams, std::marker::PhantomData);
    let info = catching_cancel(|| solver.solve(a, &pc, rhs, &mut x))?;
    let mut ax = vec![T::ZERO; n];
    a.spmv(&x, &mut ax);
    let b_norm = norm2(rhs).max(1e-300);
    let r_norm = ax.iter().zip(rhs).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt();
    let relative_residual = r_norm / b_norm;
    Ok(IterativeSolveReport {
        x,
        iterations: info.iterations,
        relative_residual,
        converged: info.reason.is_converged() && relative_residual <= (opts.rtol * 10.0).max(opts.atol),
        preconditioner_nnz,
    })
}

/// Finds `nev` interior eigenpairs of a 3D $P_2$ Nédélec curl-curl pencil `K x = λ M x` near
/// `target` via Jacobi–Davidson preconditioned by [`AmsP2`] on `K + |target| M` (no sparse LU).
pub fn eigs_curl_curl_p2_jd(
    k: &Csr<C64>,
    m: &Csr<C64>,
    topo: &AmsP2Topology,
    target: C64,
    nev: usize,
    tol: f64,
    max_it: usize,
) -> Result<CEigPairs, String> {
    use crate::linalg::jd::{JdOptions, jacobi_davidson};
    let n = k.nrows();
    if n == 0 || k.ncols() != n || m.nrows() != n || m.ncols() != n {
        return Err("eigs_curl_curl_p2_jd: K and M must be square of identical dimension".into());
    }
    if nev == 0 || nev >= n {
        return Err("eigs_curl_curl_p2_jd: nev must be between 1 and n-1".into());
    }
    progress::check()?;
    let _stage = progress::stage(format!("P2 curl-curl Jacobi–Davidson + AMS ({nev} modes, {n} unknowns)"));
    let beta = target.norm().max(1.0);
    let spd = Csr::axpby(C64::ONE, k, C64::new(beta, 0.0), m);
    let ams = AmsP2::new(&spd, topo, &AmsP2Options::default()).map_err(|e| format!("AmsP2 setup failed: {e}"))?;
    let mut opts = JdOptions::new(nev, target, true);
    opts.tol = tol;
    opts.max_it = max_it;
    opts.harmonic = false;
    opts.inner_its = 25;
    opts.m_min = (nev + 6).max(8);
    opts.m_max = (opts.m_min + 16).min(n.saturating_sub(2));
    let pc = CancelPc(&ams, std::marker::PhantomData);
    let res = catching_cancel(|| jacobi_davidson(k, Some(m), &pc, None, &opts))?
        .map_err(|e| format!("Jacobi–Davidson failed: {e}"))?;
    let vectors = (0..res.values.len()).map(|j| res.vectors.col(j).to_vec()).collect();
    Ok(CEigPairs {
        values: res.values,
        vectors,
        residuals: res.residuals,
        dim: opts.m_max,
        converged: res.converged,
    })
}

/// Matrix-free / low-memory [`ShiftedSolver`] for 3D $P_2$ Nédélec shifted curl-curl systems
/// `(K − σ M)⁻¹` solved iteratively via [`AmsP2`] + FGMRES (usable inside `ShiftInvertOperator`
/// and `krylov_schur` without allocating sparse LU factors).
pub struct IterativeShiftInvertP2<T: FaerScalar> {
    k: Csr<T>,
    m: Csr<T>,
    topo: AmsP2Topology,
    opts: IterativeCurlCurlOptions,
    shifted: Option<Csr<T>>,
    ams: Option<AmsP2<T>>,
    stats: FactorStats,
    solves: std::cell::Cell<usize>,
}

impl<T: FaerScalar> IterativeShiftInvertP2<T> {
    pub fn new(k: &Csr<T>, m: &Csr<T>, topo: AmsP2Topology, opts: IterativeCurlCurlOptions) -> Self {
        Self {
            k: k.clone(),
            m: m.clone(),
            topo,
            opts,
            shifted: None,
            ams: None,
            stats: FactorStats {
                analyse_secs: 0.0,
                symbolic_count: 0,
                numeric_secs_total: 0.0,
                numeric_secs_last: 0.0,
                numeric_count: 0,
                assemble_secs_last: 0.0,
                used: crate::linalg::direct::Used::None,
                backend: FactorBackend::None,
                cholesky_fallbacks: 0,
                solves: 0,
                factor_nnz: 0,
                solve_secs_total: 0.0,
            },
            solves: std::cell::Cell::new(0),
        }
    }
}

impl<T: FaerScalar> ShiftedSolver<T> for IterativeShiftInvertP2<T> {
    fn dim(&self) -> usize {
        self.k.nrows()
    }

    fn set_shift(&mut self, sigma: T) -> Result<(), crate::linalg::direct::SolverError> {
        let t0 = std::time::Instant::now();
        let shifted = Csr::axpby(T::ONE, &self.k, -sigma, &self.m);
        let beta = T::from_f64(sigma.modulus().max(1.0));
        let pc_mat = Csr::axpby(T::ONE, &self.k, beta, &self.m);
        let ams = AmsP2::new(&pc_mat, &self.topo, &self.opts.ams)?;
        let secs = t0.elapsed().as_secs_f64();
        self.stats.factor_nnz = ams.nnz_storage();
        self.stats.numeric_secs_last = secs;
        self.stats.numeric_secs_total += secs;
        self.stats.numeric_count += 1;
        self.shifted = Some(shifted);
        self.ams = Some(ams);
        Ok(())
    }

    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) {
        use crate::linalg::ksp::{Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances};
        let a = self.shifted.as_ref().expect("set_shift must succeed before solving");
        let ams = self.ams.as_ref().expect("set_shift must succeed before solving");
        let n = self.dim();
        let mut gmres = Gmres::<T>::new(GmresOptions {
            tol: Tolerances { rtol: self.opts.rtol, atol: self.opts.atol, max_it: self.opts.max_it, ..Default::default() },
            restart: self.opts.restart.max(2),
            side: PcSide::Right,
            orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
            flexible: self.opts.flexible,
            ..Default::default()
        });
        let mut x = vec![T::ZERO; n];
        for c in 0..nrhs {
            let col = &mut rhs[c * n..(c + 1) * n];
            let b = col.to_vec();
            let _ = gmres.solve(a, ams, &b, &mut x);
            col.copy_from_slice(&x);
        }
        self.solves.set(self.solves.get() + nrhs);
    }

    fn stats(&self) -> FactorStats {
        FactorStats { solves: self.solves.get(), ..self.stats }
    }
}

/// Fill-reducing ordering of the Cholesky, LBL^H, and custom-ordered LU paths.
/// `NYX_ORDERING=amd|auto|nd` overrides the choice.
fn fill_ordering(_n: usize) -> crate::linalg::ordering::FillOrdering {
    use crate::linalg::ordering::FillOrdering;
    match std::env::var("NYX_ORDERING").ok().as_deref() {
        Some("amd") => FillOrdering::Amd,
        Some("auto") => FillOrdering::Auto,
        Some("nd") => FillOrdering::NestedDissection,
        _ => FillOrdering::Auto,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tridiag_c64(n: usize, diag_val: C64, off_val: C64) -> CCsr {
        let mut triplets = Vec::new();
        for i in 0..n {
            triplets.push((i, i, diag_val + C64::new(i as f64 * 0.1, 0.0)));
            if i + 1 < n {
                triplets.push((i, i + 1, off_val));
                triplets.push((i + 1, i, off_val));
            }
        }
        CCsr::from_triplets(n, n, &triplets).unwrap()
    }

    #[test]
    fn direct_solver_update_matrix_reuses_symbolic_factorization() {
        let n = 6;
        let a1 = tridiag_c64(n, C64::new(4.0, 0.2), C64::new(-1.0, 0.05));
        let mut solver = DirectSolver::new(&a1, false).expect("initial factorisation");
        assert_eq!(solver.stats().symbolic_count, 1);
        assert_eq!(solver.stats().numeric_count, 1);

        let b: Vec<C64> = (0..n).map(|i| C64::new(1.0 + i as f64, -0.5 * i as f64)).collect();
        let x1 = solver.solve(&b);
        let r1 = a1.matvec(&x1);
        for i in 0..n {
            assert!((r1[i] - b[i]).norm() < 1e-12, "step 1 residual at {i}");
        }

        let a2 = tridiag_c64(n, C64::new(7.5, -0.4), C64::new(-1.2, 0.1));
        solver.update_matrix(&a2).expect("first update_matrix");
        assert_eq!(solver.stats().symbolic_count, 1, "symbolic analysis must be reused");
        assert_eq!(solver.stats().numeric_count, 2);
        let x2 = solver.solve(&b);
        let r2 = a2.matvec(&x2);
        for i in 0..n {
            assert!((r2[i] - b[i]).norm() < 1e-12, "step 2 residual at {i}");
        }

        let a3 = tridiag_c64(n, C64::new(11.0, 0.9), C64::new(-0.8, -0.2));
        solver.update_matrix(&a3).expect("second update_matrix");
        assert_eq!(solver.stats().symbolic_count, 1, "symbolic analysis must still be 1");
        assert_eq!(solver.stats().numeric_count, 3);
        let x3 = solver.solve(&b);
        let r3 = a3.matvec(&x3);
        for i in 0..n {
            assert!((r3[i] - b[i]).norm() < 1e-12, "step 3 residual at {i}");
        }
    }

    #[test]
    fn thin_svd_c64_produces_orthonormal_left_singular_vectors() {
        let col0 = vec![
            C64::new(3.0, 0.0),
            C64::new(0.0, 4.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
        ];
        let col1 = vec![
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 2.0),
            C64::new(0.0, 0.0),
        ];
        let (s, u) = thin_svd_c64(&[col0, col1]).expect("thin_svd_c64");
        assert_eq!(s.len(), 2);
        assert_eq!(u.len(), 2);
        assert!((s[0] - 5.0).abs() < 1e-12, "s[0] = {}", s[0]);
        assert!((s[1] - 2.0).abs() < 1e-12, "s[1] = {}", s[1]);

        for a in 0..2 {
            for b in 0..2 {
                let dot: C64 = u[a].iter().zip(&u[b]).map(|(x, y)| x.conj() * y).sum();
                let target = if a == b { 1.0 } else { 0.0 };
                assert!((dot - C64::new(target, 0.0)).norm() < 1e-12, "u_{a}^H u_{b} = {dot}");
            }
        }
    }

    #[test]
    fn direct_solver_stats_reports_factor_backend() {
        let a_sym = tridiag_c64(20, C64::new(4.0, 0.2), C64::new(-1.0, 0.05));
        let s_sym = DirectSolver::new(&a_sym, false).unwrap();
        assert!(
            matches!(s_sym.stats().backend, FactorBackend::LbltAmd | FactorBackend::LbltNd),
            "expected LBLT backend, got {:?}",
            s_sym.stats().backend
        );

        let mut trips = Vec::new();
        for i in 0..20 {
            trips.push((i, i, C64::new(4.0, 0.1)));
            if i + 1 < 20 {
                trips.push((i, i + 1, C64::new(-1.0, 0.2)));
                trips.push((i + 1, i, C64::new(-0.5, -0.3)));
            }
        }
        let a_unsym = CCsr::from_triplets(20, 20, &trips).unwrap();
        let s_unsym = DirectSolver::new(&a_unsym, false).unwrap();
        assert_eq!(s_unsym.stats().backend, FactorBackend::LuColamd);
    }

    #[test]
    fn small_dense_eigh_and_gen_eigh_are_bit_identical_across_pools() {
        let n = 48;
        let mut k = vec![0.0f64; n * n];
        let mut m = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = ((i * 31 + j * 17 + 7) % 97) as f64 / 97.0 - 0.5;
                let kv = if i == j { 6.0 + v } else { 0.1 * v };
                let mv = if i == j { 3.0 + 0.1 * v.abs() } else { 0.02 * v };
                k[i * n + j] = kv;
                k[j * n + i] = kv;
                m[i * n + j] = mv;
                m[j * n + i] = mv;
            }
        }
        let run_in = |nt: usize| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(nt).build().unwrap();
            pool.install(|| {
                let (w1, v1) = dense_eigh(&k, n).unwrap();
                let (w2, v2) = dense_gen_eigh(&k, &m, n).unwrap();
                (w1, v1, w2, v2)
            })
        };
        let (w1_1, v1_1, w2_1, v2_1) = run_in(1);
        let (w1_8, v1_8, w2_8, v2_8) = run_in(8);
        assert_eq!(w1_1.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), w1_8.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
        assert_eq!(w2_1.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), w2_8.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
        for (a, b) in v1_1.iter().flatten().zip(v1_8.iter().flatten()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        for (a, b) in v2_1.iter().flatten().zip(v2_8.iter().flatten()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}
