//! Nonlinear solvers for `F(x) = 0` (port item L1; the role of PETSc `SNES`): Newton with a
//! backtracking line search and an optional Jacobian lag, Jacobian-free Newton–Krylov (JFNK) with
//! Eisenstat–Walker forcing, Picard iteration and Anderson acceleration of fixed-point maps.
//!
//! Uses: temperature-dependent conductivity, Kerr and thermo-optic nonlinearities, photothermal
//! coupling. All solvers are generic over real and complex scalars; for complex `x`, `F` must be
//! *complex-analytic* (its Jacobian a complex matrix); non-analytic maps such as `|u|²u` are solved as a
//! real system of twice the size.
//!
//! * [`newton`]: `J(x_k) d = −F(x_k)` by sparse direct factorisation (one symbolic analysis for the
//!   whole solve, values replaced every Jacobian evaluation), step `x ← x + α d` with `α` from quadratic
//!   backtracking on `½‖F‖²` (Armijo, `c = 1e-4`) or `α = 1`. `jacobian_lag = L` reuses the Jacobian and
//!   its factorisation for `L` iterations (modified Newton: fewer factorisations, linear convergence).
//! * [`jfnk`]: the Jacobian-vector product is a forward difference `(F(x + εv) − F(x))/ε`,
//!   the linear system is solved by GMRES to the Eisenstat–Walker tolerance, and any
//!   [`Preconditioner`] (an LU or AMG of an approximate Jacobian) can be supplied.
//! * [`picard`], [`anderson`]: fixed-point maps `x ← G(x)`; Anderson (type II, Walker & Ni 2011) mixes
//!   the last `m` residuals by a least-squares problem solved with QR (modified Gram–Schmidt).
//!
//! Convergence tests follow PETSc: `‖F‖ ≤ max(rtol·‖F₀‖, atol)`; reasons are reported, never hidden.

use crate::solvers::csr::Csr;
use crate::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};
use crate::solvers::faer_glue::FaerScalar;
use crate::solvers::ksp::{ConvergedReason, Gmres, GmresOptions, PcSide, Tolerances};
use crate::solvers::linop::{FnOperator, Identity, Preconditioner};
use crate::solvers::scalar::Scalar;
use crate::solvers::vecops::norm2;
use std::cell::RefCell;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reason {
    ConvergedRtol,
    ConvergedAtol,
    /// The step became smaller than `stol·‖x‖` *and* `‖F‖ ≤ √rtol·‖F₀‖` (unlike PETSc's `snes_stol` test, a negligible step
    /// with a large residual is reported as [`Reason::DivergedStagnation`]).
    ConvergedStep,
    DivergedMaxIt,
    /// The line search could not reduce `‖F‖`.
    DivergedLineSearch,
    /// The step became negligible while `‖F‖` is still far from the tolerance (PETSc would call this converged).
    DivergedStagnation,
    /// `F` returned NaN or infinity.
    DivergedFunctionNorm,
    /// The linear solve failed (singular Jacobian or Krylov breakdown).
    DivergedLinearSolve,
}

impl Reason {
    pub fn is_converged(self) -> bool {
        matches!(self, Reason::ConvergedRtol | Reason::ConvergedAtol | Reason::ConvergedStep)
    }
}

#[derive(Debug, Clone)]
pub struct NonlinearResult<T> {
    pub x: Vec<T>,
    pub reason: Reason,
    pub iterations: usize,
    /// `‖F(x_k)‖` for `k = 0 ..= iterations`.
    pub residual_norms: Vec<f64>,
    pub function_evals: usize,
    pub jacobian_evals: usize,
    /// Numeric factorisations (Newton) and symbolic analyses.
    pub factorizations: usize,
    pub symbolic_analyses: usize,
    /// Krylov iterations in total (JFNK).
    pub linear_iterations: usize,
    /// Backtracking steps taken in total.
    pub backtracks: usize,
}

impl<T> NonlinearResult<T> {
    pub fn converged(&self) -> bool {
        self.reason.is_converged()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineSearch {
    /// Full steps.
    None,
    /// Quadratic backtracking with Armijo constant `c` and a floor on the step length.
    Backtracking { c: f64, alpha_min: f64 },
}

impl Default for LineSearch {
    fn default() -> Self {
        LineSearch::Backtracking { c: 1e-4, alpha_min: 1e-10 }
    }
}

#[derive(Debug, Clone)]
pub struct NewtonOptions {
    pub rtol: f64,
    pub atol: f64,
    /// Stop when `‖αd‖ ≤ stol·‖x‖` (default 1e-12, a stagnation guard).
    pub stol: f64,
    pub max_it: usize,
    pub line_search: LineSearch,
    /// Jacobian re-evaluated and refactorised every `jacobian_lag` iterations (1 = every iteration).
    pub jacobian_lag: usize,
    pub definiteness: Definiteness,
}

impl Default for NewtonOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-8,
            atol: 1e-50,
            stol: 1e-12,
            max_it: 50,
            line_search: LineSearch::default(),
            jacobian_lag: 1,
            definiteness: Definiteness::General,
        }
    }
}

fn axpy_into<T: Scalar>(x: &[T], alpha: f64, d: &[T], out: &mut [T]) {
    for ((o, &xi), &di) in out.iter_mut().zip(x).zip(d) {
        *o = xi + di.scale(alpha);
    }
}

struct Tracker<T> {
    norms: Vec<f64>,
    fevals: usize,
    _p: std::marker::PhantomData<T>,
}

fn converged(norm: f64, f0: f64, rtol: f64, atol: f64) -> Option<Reason> {
    if norm <= atol {
        Some(Reason::ConvergedAtol)
    } else if norm <= rtol * f0 {
        Some(Reason::ConvergedRtol)
    } else {
        None
    }
}

/// Quadratic backtracking on `φ(α) = ½‖F(x + αd)‖²`; returns `(α, ‖F‖, backtracks)` or `None` on failure.
#[allow(clippy::too_many_arguments)]
fn backtrack<T: Scalar>(
    f: &mut dyn FnMut(&[T], &mut [T]),
    x: &[T],
    d: &[T],
    fnorm: f64,
    ls: LineSearch,
    trial: &mut [T],
    ftrial: &mut [T],
    fevals: &mut usize,
) -> Option<(f64, f64, usize)> {
    match ls {
        LineSearch::None => {
            axpy_into(x, 1.0, d, trial);
            f(trial, ftrial);
            *fevals += 1;
            let nn = norm2(ftrial);
            if nn.is_finite() { Some((1.0, nn, 0)) } else { None }
        }
        LineSearch::Backtracking { c, alpha_min } => {
            let phi0 = fnorm * fnorm;
            let (mut alpha, mut bt) = (1.0f64, 0usize);
            loop {
                axpy_into(x, alpha, d, trial);
                f(trial, ftrial);
                *fevals += 1;
                let nn = norm2(ftrial);
                // exact Newton direction: ‖F + αJd‖ = (1 − α)‖F‖, so ½-merit slope is −‖F‖²
                if nn.is_finite() && nn * nn <= (1.0 - 2.0 * c * alpha) * phi0 {
                    return Some((alpha, nn, bt));
                }
                if alpha <= alpha_min {
                    return None;
                }
                // minimiser of the quadratic through φ(0), φ'(0) = −2φ0 (in ‖F‖² units) and φ(α)
                let phi_a = if nn.is_finite() { nn * nn } else { f64::INFINITY };
                let denom = phi_a - phi0 + 2.0 * phi0 * alpha;
                let next = if denom > 0.0 && phi_a.is_finite() { phi0 * alpha * alpha / denom } else { 0.1 * alpha };
                alpha = next.clamp(0.1 * alpha, 0.5 * alpha).max(alpha_min * 0.5);
                bt += 1;
            }
        }
    }
}

/// Newton's method with a sparse direct linear solve (see the module docs).
pub fn newton<T: FaerScalar>(
    f: &mut dyn FnMut(&[T], &mut [T]),
    jac: &mut dyn FnMut(&[T]) -> Csr<T>,
    x0: &[T],
    opts: &NewtonOptions,
) -> Result<NonlinearResult<T>, SolverError> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let (mut fx, mut trial, mut ftrial, mut d) =
        (vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n]);
    let mut tr = Tracker::<T> { norms: Vec::new(), fevals: 0, _p: Default::default() };
    f(&x, &mut fx);
    tr.fevals += 1;
    let mut fnorm = norm2(&fx);
    let f0 = fnorm;
    tr.norms.push(fnorm);
    let (mut jevals, mut backtracks, mut iterations) = (0usize, 0usize, 0usize);
    let mut fac: Option<ShiftedFactorization<T>> = None;
    let mut since_jac = usize::MAX;
    let mut reason = if !fnorm.is_finite() {
        Some(Reason::DivergedFunctionNorm)
    } else {
        converged(fnorm, f0, opts.rtol, opts.atol)
    };
    while reason.is_none() {
        if iterations >= opts.max_it {
            reason = Some(Reason::DivergedMaxIt);
            break;
        }
        if since_jac >= opts.jacobian_lag.max(1) {
            let j = jac(&x);
            jevals += 1;
            let ff = match fac.as_mut() {
                Some(ff) => {
                    ff.replace_matrices(&j, None)?;
                    ff
                }
                None => fac.insert(ShiftedFactorization::new(&j, None, opts.definiteness)?),
            };
            if ff.set_shift(T::ZERO).is_err() {
                reason = Some(Reason::DivergedLinearSolve);
                break;
            }
            since_jac = 0;
        }
        since_jac += 1;
        for (di, fi) in d.iter_mut().zip(&fx) {
            *di = -*fi;
        }
        fac.as_ref().unwrap().solve_block_in_place(&mut d, 1);
        if d.iter().any(|v| !v.is_finite()) {
            reason = Some(Reason::DivergedLinearSolve);
            break;
        }
        let Some((alpha, new_norm, bt)) =
            backtrack(f, &x, &d, fnorm, opts.line_search, &mut trial, &mut ftrial, &mut tr.fevals)
        else {
            reason = Some(if ftrial.iter().any(|v| !v.is_finite()) {
                Reason::DivergedFunctionNorm
            } else {
                Reason::DivergedLineSearch
            });
            break;
        };
        backtracks += bt;
        let step = alpha * norm2(&d);
        x.copy_from_slice(&trial);
        fx.copy_from_slice(&ftrial);
        fnorm = new_norm;
        iterations += 1;
        tr.norms.push(fnorm);
        reason = converged(fnorm, f0, opts.rtol, opts.atol);
        if reason.is_none() && step <= opts.stol * norm2(&x).max(f64::MIN_POSITIVE) {
            reason =
                Some(if fnorm <= opts.rtol.sqrt() * f0 { Reason::ConvergedStep } else { Reason::DivergedStagnation });
        }
    }
    let (fz, sy) = fac.as_ref().map_or((0, 0), |ff| (ff.stats().numeric_count, ff.stats().symbolic_count));
    Ok(NonlinearResult {
        x,
        reason: reason.unwrap(),
        iterations,
        residual_norms: tr.norms,
        function_evals: tr.fevals,
        jacobian_evals: jevals,
        factorizations: fz,
        symbolic_analyses: sy,
        linear_iterations: 0,
        backtracks,
    })
}

// ---------------------------------------------------------------------------------------------
// JFNK
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Forcing {
    Fixed(f64),
    /// Eisenstat–Walker choice 2: `η_k = γ (‖F_k‖/‖F_{k−1}‖)^α`, `γ = 1`, `α = (1+√5)/2`, clipped to `[1e-6, 0.9]`
    /// and safeguarded against a too fast decrease.
    EisenstatWalker2,
}

#[derive(Debug, Clone)]
pub struct JfnkOptions {
    pub rtol: f64,
    pub atol: f64,
    pub max_it: usize,
    pub line_search: LineSearch,
    pub forcing: Forcing,
    /// GMRES restart and iteration cap per Newton step.
    pub restart: usize,
    pub max_linear_it: usize,
}

impl Default for JfnkOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-8,
            atol: 1e-50,
            max_it: 50,
            line_search: LineSearch::default(),
            forcing: Forcing::EisenstatWalker2,
            restart: 30,
            max_linear_it: 100,
        }
    }
}

/// Jacobian-free Newton–Krylov. `pc` approximates `J⁻¹` (identity if `None`), is applied from the *right* (so the
/// forcing tolerance bounds the true linear residual) and is fixed for the solve.
pub fn jfnk<T: Scalar>(
    f: &mut dyn FnMut(&[T], &mut [T]),
    x0: &[T],
    pc: Option<&dyn Preconditioner<T>>,
    opts: &JfnkOptions,
) -> NonlinearResult<T> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let (mut fx, mut trial, mut ftrial, mut d) =
        (vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n]);
    let f_cell = RefCell::new(&mut *f);
    let fevals_cell = RefCell::new(0usize);
    let mut tr = Tracker::<T> { norms: Vec::new(), fevals: 0, _p: Default::default() };
    (f_cell.borrow_mut())(&x, &mut fx);
    *fevals_cell.borrow_mut() += 1;
    let mut fnorm = norm2(&fx);
    let (f0, mut fprev) = (fnorm, fnorm);
    tr.norms.push(fnorm);
    let (mut linear_its, mut backtracks, mut iterations) = (0usize, 0usize, 0usize);
    let mut eta = match opts.forcing {
        Forcing::Fixed(e) => e,
        Forcing::EisenstatWalker2 => 0.3,
    };
    let mut reason = if !fnorm.is_finite() {
        Some(Reason::DivergedFunctionNorm)
    } else {
        converged(fnorm, f0, opts.rtol, opts.atol)
    };
    let ident = Identity;
    let pcr: &dyn Preconditioner<T> = pc.unwrap_or(&ident);
    while reason.is_none() {
        if iterations >= opts.max_it {
            reason = Some(Reason::DivergedMaxIt);
            break;
        }
        if iterations > 0 && opts.forcing == Forcing::EisenstatWalker2 {
            {
                let ratio = fnorm / fprev.max(f64::MIN_POSITIVE);
                let alpha = 0.5 * (1.0 + 5.0f64.sqrt());
                let mut e = ratio.powf(alpha);
                // safeguard: do not let the tolerance drop faster than the previous one allows
                let prev_pow = eta.powf(alpha);
                if prev_pow > 0.1 {
                    e = e.max(prev_pow);
                }
                eta = e.clamp(1e-6, 0.9);
            }
        }
        let (xs, fxs) = (x.clone(), fx.clone());
        let xnorm = norm2(&xs);
        // J v ≈ (F(x + εv) − F(x)) / ε
        let jv = |v: &[T], out: &mut [T]| {
            let vn = norm2(v);
            if vn == 0.0 {
                out.fill(T::ZERO);
                return;
            }
            let eps = f64::EPSILON.sqrt() * (1.0 + xnorm) / vn;
            let pert: Vec<T> = xs.iter().zip(v).map(|(&a, &b)| a + b.scale(eps)).collect();
            (f_cell.borrow_mut())(&pert, out);
            *fevals_cell.borrow_mut() += 1;
            for (o, &f0i) in out.iter_mut().zip(&fxs) {
                *o = (*o - f0i).scale(1.0 / eps);
            }
        };
        let op = FnOperator::new(n, jv);
        let rhs: Vec<T> = fx.iter().map(|&v| -v).collect();
        let mut gm = Gmres::<T>::new(GmresOptions {
            tol: Tolerances { rtol: eta, atol: 1e-50, max_it: opts.max_linear_it, ..Default::default() },
            restart: opts.restart,
            // right preconditioning: the Krylov residual is the *true* linear residual, which is what the
            // forcing term `η` (inexact-Newton theory) refers to
            side: PcSide::Right,
            ..Default::default()
        });
        d.fill(T::ZERO);
        let info = gm.solve(&op, pcr, &rhs, &mut d);
        linear_its += info.iterations;
        if matches!(info.reason, ConvergedReason::DivergedBreakdown | ConvergedReason::DivergedNanOrInf)
            || d.iter().any(|v| !v.is_finite())
        {
            reason = Some(Reason::DivergedLinearSolve);
            break;
        }
        let mut fe = 0usize;
        let ls = {
            let mut fb = |a: &[T], b: &mut [T]| (f_cell.borrow_mut())(a, b);
            backtrack(&mut fb, &x, &d, fnorm, opts.line_search, &mut trial, &mut ftrial, &mut fe)
        };
        *fevals_cell.borrow_mut() += fe;
        let Some((alpha, new_norm, bt)) = ls else {
            reason = Some(if ftrial.iter().any(|v| !v.is_finite()) {
                Reason::DivergedFunctionNorm
            } else {
                Reason::DivergedLineSearch
            });
            break;
        };
        let _ = alpha;
        backtracks += bt;
        x.copy_from_slice(&trial);
        fx.copy_from_slice(&ftrial);
        fprev = fnorm;
        fnorm = new_norm;
        iterations += 1;
        tr.norms.push(fnorm);
        reason = converged(fnorm, f0, opts.rtol, opts.atol);
    }
    tr.fevals = *fevals_cell.borrow();
    NonlinearResult {
        x,
        reason: reason.unwrap(),
        iterations,
        residual_norms: tr.norms,
        function_evals: tr.fevals,
        jacobian_evals: 0,
        factorizations: 0,
        symbolic_analyses: 0,
        linear_iterations: linear_its,
        backtracks,
    }
}

// ---------------------------------------------------------------------------------------------
// Fixed-point iterations
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct FixedPointOptions {
    pub rtol: f64,
    pub atol: f64,
    pub max_it: usize,
    /// Relaxation `β` (1 = plain Picard step `x ← G(x)`).
    pub beta: f64,
    /// Anderson depth `m` (0 = plain Picard).
    pub depth: usize,
}

impl Default for FixedPointOptions {
    fn default() -> Self {
        Self { rtol: 1e-8, atol: 1e-50, max_it: 500, beta: 1.0, depth: 0 }
    }
}

/// Picard iteration `x ← x + β (G(x) − x)`.
pub fn picard<T: Scalar>(g: &mut dyn FnMut(&[T], &mut [T]), x0: &[T], opts: &FixedPointOptions) -> NonlinearResult<T> {
    anderson(g, x0, &FixedPointOptions { depth: 0, ..opts.clone() })
}

/// Anderson acceleration of `x ← G(x)` (type II, depth `opts.depth`); depth 0 is Picard. The residual
/// is `f_k = G(x_k) − x_k` and the test is `‖f_k‖ ≤ max(rtol‖f_0‖, atol)`.
pub fn anderson<T: Scalar>(
    g: &mut dyn FnMut(&[T], &mut [T]),
    x0: &[T],
    opts: &FixedPointOptions,
) -> NonlinearResult<T> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut gx = vec![T::ZERO; n];
    let mut fk = vec![T::ZERO; n];
    let mut norms = Vec::new();
    let (mut dx, mut df): (Vec<Vec<T>>, Vec<Vec<T>>) = (Vec::new(), Vec::new());
    let (mut x_prev, mut f_prev): (Option<Vec<T>>, Option<Vec<T>>) = (None, None);
    let mut fevals = 0usize;
    let mut iterations = 0usize;
    let mut f0 = 0.0;
    let reason = loop {
        g(&x, &mut gx);
        fevals += 1;
        for i in 0..n {
            fk[i] = gx[i] - x[i];
        }
        let fnorm = norm2(&fk);
        if iterations == 0 {
            f0 = fnorm;
        }
        norms.push(fnorm);
        if !fnorm.is_finite() {
            break Reason::DivergedFunctionNorm;
        }
        if let Some(r) = converged(fnorm, f0, opts.rtol, opts.atol) {
            break r;
        }
        if iterations >= opts.max_it {
            break Reason::DivergedMaxIt;
        }
        // update the difference histories
        if let (Some(xp), Some(fp)) = (&x_prev, &f_prev) {
            dx.push(x.iter().zip(xp).map(|(&a, &b)| a - b).collect());
            df.push(fk.iter().zip(fp).map(|(&a, &b)| a - b).collect());
            if dx.len() > opts.depth {
                dx.remove(0);
                df.remove(0);
            }
        }
        x_prev = Some(x.clone());
        f_prev = Some(fk.clone());
        let mut next: Vec<T> = (0..n).map(|i| x[i] + fk[i].scale(opts.beta)).collect();
        if opts.depth > 0 && !df.is_empty() {
            // least squares  min_γ ‖f_k − ΔF γ‖  by modified Gram–Schmidt QR, dropping dependent columns
            let m = df.len();
            let mut q: Vec<Vec<T>> = Vec::with_capacity(m);
            let mut r = vec![vec![T::ZERO; m]; m];
            let mut keep: Vec<usize> = Vec::new();
            for (j, col) in df.iter().enumerate() {
                let mut v = col.clone();
                let scale = norm2(&v).max(f64::MIN_POSITIVE);
                for (qi, &ki) in q.iter().zip(&keep) {
                    let h = crate::solvers::vecops::dotc(qi, &v);
                    r[ki][j] = h;
                    for (vv, qq) in v.iter_mut().zip(qi) {
                        *vv -= h * *qq;
                    }
                }
                let nv = norm2(&v);
                if nv > 1e-10 * scale {
                    v.iter_mut().for_each(|e| *e = e.scale(1.0 / nv));
                    r[j][j] = T::from_f64(nv);
                    q.push(v);
                    keep.push(j);
                }
            }
            let rhs: Vec<T> = q.iter().map(|qi| crate::solvers::vecops::dotc(qi, &fk)).collect();
            // back substitution on the kept columns
            let k = keep.len();
            let mut gamma = vec![T::ZERO; m];
            for a in (0..k).rev() {
                let ia = keep[a];
                let mut s = rhs[a];
                for b in a + 1..k {
                    s -= r[ia][keep[b]] * gamma[keep[b]];
                }
                gamma[ia] = s / r[ia][ia];
            }
            for (j, gj) in gamma.iter().enumerate() {
                for i in 0..n {
                    next[i] -= (dx[j][i] + df[j][i].scale(opts.beta)) * *gj;
                }
            }
        }
        x = next;
        iterations += 1;
    };
    NonlinearResult {
        x,
        reason,
        iterations,
        residual_norms: norms,
        function_evals: fevals,
        jacobian_evals: 0,
        factorizations: 0,
        symbolic_analyses: 0,
        linear_iterations: 0,
        backtracks: 0,
    }
}
