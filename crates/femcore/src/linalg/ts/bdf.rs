//! Variable-step, variable-order BDF (orders 1–5) for `M y' + K y = f(t)`.
//!
//! The formulation is the *quasi-constant step size* one of Shampine & Reichelt (1997, "The
//! MATLAB ODE suite"): the solution history is kept as backward differences `D_k = ∇^k y_n`
//! of an equally spaced grid with the current step `h`; changing `h` is a linear transformation of
//! `D` (matrices `R`, `U` of the paper), so the formulas keep fixed coefficients. Written from the
//! paper, not from any library's code.
//!
//! One step of order `p`, with `γ_k = Σ_{j≤k} 1/j`:
//!
//! ```text
//! y_pred = Σ_{k=0..p} D_k              (extrapolation)
//! ψ      = Σ_{k=1..p} γ_k D_k / γ_p,   c = h / γ_p
//! (K + M/c) d = f(t_{n+1}) − K y_pred − M ψ / c,       y_{n+1} = y_pred + d
//! err    = d / (p + 1)                  (local error estimate)
//! ```
//!
//! (a *linear* problem, so the "Newton" iteration is a single solve with the exact matrix). The
//! factorisation of `K + M/c` is reused for as long as `c` is unchanged; `c` changes only when the
//! step size or the order changes, which the controller allows at most once per `order + 1` steps.
//!
//! Step/order selection follows the paper: after `order + 1` steps at the same `h` the error is
//! estimated for orders `p − 1`, `p`, `p + 1` and the best predicted step is taken (growth capped at
//! 10×, shrink floor 0.2×, safety 0.9); a rejected step shrinks `h` and keeps the order.
//!
//! Output at arbitrary times uses the interpolating polynomial of the history (degree = order).

use super::{Solution, TsError, TsStats, scaled_rms};
use crate::linalg::csr::Csr;
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::linalg::faer_glue::FaerScalar;
use std::time::Instant;

const MAX_ORDER: usize = 5;
const SAFETY: f64 = 0.9;
const MIN_FACTOR: f64 = 0.2;
const MAX_FACTOR: f64 = 10.0;

#[derive(Debug, Clone)]
pub struct BdfOptions {
    pub rtol: f64,
    pub atol: f64,
    /// Highest order used (1..=5).
    pub max_order: usize,
    /// Initial step (otherwise chosen from `y'(t0)`).
    pub first_step: Option<f64>,
    pub max_step: f64,
    /// Smallest step before giving up, relative to `max(1, |t_end − t0|)`.
    pub min_step_rel: f64,
    /// Constant step size, no error control (the order still ramps up one per `order + 1` steps,
    /// up to `max_order`). For convergence studies.
    pub fixed_step: Option<f64>,
    pub max_steps: usize,
    /// `K + M/c` is SPD? (`PositiveDefinite` tries Cholesky first.)
    pub definiteness: Definiteness,
    /// Reuse the last factorisation of `K + M/c₀` for `c` within `±reuse_tolerance·c₀` by iterative
    /// refinement (each refinement costs one solve and two sparse products; it converges with ratio
    /// `≤ |c₀/c − 1|` for SPD `K`), refactorising only outside that band or when refinement does not
    /// converge in `MAX_REFINE` solves. `0.0` (default) = always factorise exactly at the current `c`:
    /// measured on 3-D heat (64 000 unknowns) a band of ±30 % saved 9 of 26 factorisations but lost
    /// 8 % of the wall time to the refinement solves; it only pays when a factorisation costs far
    /// more than ~30 solves.
    pub reuse_tolerance: f64,
}

impl BdfOptions {
    pub fn new(rtol: f64, atol: f64) -> Self {
        Self {
            rtol,
            atol,
            max_order: MAX_ORDER,
            first_step: None,
            max_step: f64::INFINITY,
            min_step_rel: 1e-14,
            fixed_step: None,
            max_steps: 1_000_000,
            definiteness: Definiteness::General,
            reuse_tolerance: 0.0,
        }
    }
}

/// `R` of the paper: `R[i][j] = Π_{r=1..i} (r − 1 − factor·j)/r` (row 0 is all ones), `0 ≤ i, j ≤ order`.
fn r_matrix(order: usize, factor: f64) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; order + 1]; order + 1];
    for j in 0..=order {
        m[0][j] = 1.0;
    }
    for i in 1..=order {
        for j in 1..=order {
            m[i][j] = m[i - 1][j] * ((i as f64 - 1.0 - factor * j as f64) / i as f64);
        }
        m[i][0] = 0.0; // column 0 is the unit vector (the value y itself does not change)
    }
    m[0][0] = 1.0;
    m
}

/// Rescales the backward differences `D[0..=order]` to the step `factor · h`.
fn change_step<T: FaerScalar>(d: &mut [Vec<T>], order: usize, factor: f64) {
    let r = r_matrix(order, factor);
    let u = r_matrix(order, 1.0);
    // RU = R·U ;  D_new[i] = Σ_j RU[j][i] D_old[j]
    let mut ru = vec![vec![0.0; order + 1]; order + 1];
    for i in 0..=order {
        for j in 0..=order {
            ru[i][j] = (0..=order).map(|k| r[i][k] * u[k][j]).sum();
        }
    }
    let n = d[0].len();
    let old: Vec<Vec<T>> = d[..=order].to_vec();
    for i in 0..=order {
        for e in 0..n {
            let mut s = T::ZERO;
            for j in 0..=order {
                s += old[j][e].scale(ru[j][i]);
            }
            d[i][e] = s;
        }
    }
}

/// Value of the history polynomial at `t_new + x·h` (`x ≤ 0` interpolates, small `x > 0` extrapolates).
fn history_value<T: FaerScalar>(d: &[Vec<T>], order: usize, x: f64, out: &mut [T]) {
    out.copy_from_slice(&d[0]);
    let mut coef = 1.0;
    for (k, dk) in d.iter().enumerate().take(order + 1).skip(1) {
        coef *= (x + (k - 1) as f64) / k as f64;
        for (o, v) in out.iter_mut().zip(dk) {
            *o += v.scale(coef);
        }
    }
}

struct IterationMatrix<T: FaerScalar> {
    fac: ShiftedFactorization<T>,
    /// `c` of the current numeric factorisation.
    c: Option<f64>,
    reuse_tolerance: f64,
}

const MAX_REFINE: usize = 4;

impl<T: FaerScalar> IterationMatrix<T> {
    fn refactor(&mut self, c: f64) -> Result<(), TsError> {
        // K + M/c = K − σ M with σ = −1/c
        self.fac.set_shift(T::from_f64(-1.0 / c))?;
        self.c = Some(c);
        Ok(())
    }

    /// Solves `(K + M/c) d = rhs` into `d`; `scale_ref`/`rtol`/`atol` define the refinement target
    /// (a correction below 1 % of the integration tolerance).
    #[allow(clippy::too_many_arguments)]
    fn solve(
        &mut self,
        k: &Csr<T>,
        m: &Csr<T>,
        c: f64,
        rhs: &[T],
        d: &mut [T],
        scale_ref: &[T],
        rtol: f64,
        atol: f64,
    ) -> Result<(), TsError> {
        let n = rhs.len();
        let stale = match self.c {
            Some(c0) => (c0 - c).abs() > 1e-13 * c.abs() && (c0 - c).abs() <= self.reuse_tolerance * c0,
            None => false,
        };
        if !stale {
            if self.c.is_none_or(|c0| (c0 - c).abs() > 1e-13 * c.abs()) {
                self.refactor(c)?;
            }
            d.copy_from_slice(rhs);
            self.fac.solve_block_in_place(d, 1)?;
            return Ok(());
        }
        // modified fixed-point iteration with the stale factorisation
        d.copy_from_slice(rhs);
        self.fac.solve_block_in_place(d, 1)?;
        let (mut res, mut tmp) = (vec![T::ZERO; n], vec![T::ZERO; n]);
        for _ in 0..MAX_REFINE {
            k.spmv(d, &mut res);
            m.spmv(d, &mut tmp);
            for i in 0..n {
                res[i] = rhs[i] - res[i] - tmp[i].scale(1.0 / c);
            }
            self.fac.solve_block_in_place(&mut res, 1)?;
            // the stale factor solves (K + M/c₀): the correction is (K + M/c₀)⁻¹ r
            let size = scaled_rms(&res, scale_ref, rtol, atol);
            for i in 0..n {
                d[i] += res[i];
            }
            if size < 1e-2 {
                return Ok(());
            }
        }
        // did not converge: factorise at the current c and solve exactly
        self.refactor(c)?;
        d.copy_from_slice(rhs);
        self.fac.solve_block_in_place(d, 1)?;
        Ok(())
    }
}

/// Integrates `M y' + K y = f(t)` from `t0` to `t_end`, returning `y` at `t_eval` (sorted, inside
/// `[t0, t_end]`; empty = every accepted step). `force(t, f)` writes `f(t)`.
#[allow(clippy::too_many_arguments)]
pub fn bdf_solve<T: FaerScalar>(
    m: &Csr<T>,
    k: &Csr<T>,
    force: &mut dyn FnMut(f64, &mut [T]),
    t0: f64,
    t_end: f64,
    y0: &[T],
    t_eval: &[f64],
    opts: &BdfOptions,
) -> Result<Solution<T>, TsError> {
    let n = y0.len();
    if m.nrows() != n || m.ncols() != n || k.nrows() != n || k.ncols() != n {
        return Err(TsError::InvalidArgument("M and K must be n×n with n = len(y0)"));
    }
    if !(t_end > t0) {
        return Err(TsError::InvalidArgument("t_end must exceed t0"));
    }
    if opts.max_order == 0 || opts.max_order > MAX_ORDER {
        return Err(TsError::InvalidArgument("max_order must be in 1..=5"));
    }
    if t_eval.windows(2).any(|w| w[1] < w[0]) || t_eval.iter().any(|&t| t < t0 || t > t_end) {
        return Err(TsError::InvalidArgument("t_eval must be sorted and inside [t0, t_end]"));
    }
    let started = Instant::now();
    let span = (t_end - t0).abs().max(1.0);
    let min_step = opts.min_step_rel * span;
    let gamma: Vec<f64> = (0..=MAX_ORDER + 1)
        .scan(0.0, |s, j| {
            if j > 0 {
                *s += 1.0 / j as f64;
            }
            Some(*s)
        })
        .collect();

    let mut itm = IterationMatrix {
        fac: ShiftedFactorization::new(k, Some(m), opts.definiteness)?,
        c: None,
        reuse_tolerance: opts.reuse_tolerance,
    };
    let mut mfac = ShiftedFactorization::new(m, None, opts.definiteness)?;
    mfac.set_shift(T::ZERO)?;

    let mut stats = TsStats::default();
    let mut f = vec![T::ZERO; n];
    let mut tmp = vec![T::ZERO; n];

    // y'(t0) = M⁻¹ (f(t0) − K y0)
    let mut ydot0 = vec![T::ZERO; n];
    force(t0, &mut f);
    stats.forcing_evals += 1;
    k.spmv(y0, &mut tmp);
    for i in 0..n {
        ydot0[i] = f[i] - tmp[i];
    }
    mfac.solve_block_in_place(&mut ydot0, 1)?;

    // first step
    let mut h = match (opts.fixed_step, opts.first_step) {
        (Some(h), _) | (None, Some(h)) => h,
        _ => {
            let scale0: Vec<T> = y0.to_vec();
            let d0 = scaled_rms(y0, &scale0, opts.rtol, opts.atol);
            let d1 = scaled_rms(&ydot0, &scale0, opts.rtol, opts.atol);
            let h0 = if d0 < 1e-5 || d1 < 1e-5 { 1e-6 } else { 0.01 * d0 / d1 };
            // one explicit Euler probe for the second derivative
            let y1: Vec<T> = y0.iter().zip(&ydot0).map(|(&a, &b)| a + b.scale(h0)).collect();
            force(t0 + h0, &mut f);
            stats.forcing_evals += 1;
            k.spmv(&y1, &mut tmp);
            let mut ydot1: Vec<T> = f.iter().zip(&tmp).map(|(&a, &b)| a - b).collect();
            mfac.solve_block_in_place(&mut ydot1, 1)?;
            let diff: Vec<T> = ydot1.iter().zip(&ydot0).map(|(&a, &b)| a - b).collect();
            let d2 = scaled_rms(&diff, &scale0, opts.rtol, opts.atol) / h0;
            let h1 = if d1.max(d2) <= 1e-15 { (h0 * 1e-3).max(1e-6) } else { (0.01 / d1.max(d2)).sqrt() };
            (100.0 * h0).min(h1)
        }
    };
    h = h.min(opts.max_step).min(t_end - t0);
    if !(h > 0.0) {
        return Err(TsError::InvalidArgument("step size must be positive"));
    }

    let mut d: Vec<Vec<T>> = vec![vec![T::ZERO; n]; MAX_ORDER + 3];
    d[0].copy_from_slice(y0);
    for i in 0..n {
        d[1][i] = ydot0[i].scale(h);
    }
    let mut order = 1usize;
    let mut n_equal = 0usize;
    let mut t = t0;

    let mut out_t = Vec::new();
    let mut out_y: Vec<Vec<T>> = Vec::new();
    let mut next_eval = 0usize;
    // outputs at t0 itself
    while next_eval < t_eval.len() && t_eval[next_eval] <= t0 {
        out_t.push(t_eval[next_eval]);
        out_y.push(y0.to_vec());
        next_eval += 1;
    }
    if t_eval.is_empty() {
        out_t.push(t0);
        out_y.push(y0.to_vec());
    }

    let mut y_pred = vec![T::ZERO; n];
    let mut psi = vec![T::ZERO; n];
    let mut dcorr = vec![T::ZERO; n];
    let mut y_new = vec![T::ZERO; n];
    let mut err = vec![T::ZERO; n];
    let mut rhs = vec![T::ZERO; n];
    let mut steps_taken = 0usize;

    while t < t_end {
        steps_taken += 1;
        if steps_taken > opts.max_steps {
            return Err(TsError::TooManySteps(t));
        }
        if h < min_step {
            return Err(TsError::StepTooSmall(t));
        }
        // do not step beyond the end (or leave a sliver)
        if t + h > t_end || (opts.fixed_step.is_none() && t_end - (t + h) < 1e-12 * span) {
            let h_new = t_end - t;
            if (h_new - h).abs() > 1e-15 * h {
                change_step(&mut d, order, h_new / h);
                h = h_new;
                n_equal = 0;
            }
        }
        let t_new = t + h;

        // predictor and history combination
        y_pred.copy_from_slice(&d[0]);
        for dk in d.iter().take(order + 1).skip(1) {
            for (p, v) in y_pred.iter_mut().zip(dk) {
                *p += *v;
            }
        }
        for e in 0..n {
            let mut s = T::ZERO;
            for kk in 1..=order {
                s += d[kk][e].scale(gamma[kk]);
            }
            psi[e] = s.scale(1.0 / gamma[order]);
        }
        let c = h / gamma[order];
        force(t_new, &mut f);
        stats.forcing_evals += 1;
        k.spmv(&y_pred, &mut tmp);
        m.spmv(&psi, &mut rhs);
        for e in 0..n {
            rhs[e] = f[e] - tmp[e] - rhs[e].scale(1.0 / c);
        }
        itm.solve(k, m, c, &rhs, &mut dcorr, &y_pred, opts.rtol, opts.atol)?;
        for e in 0..n {
            y_new[e] = y_pred[e] + dcorr[e];
            err[e] = dcorr[e].scale(1.0 / (order as f64 + 1.0));
        }
        let err_norm = scaled_rms(&err, &y_new, opts.rtol, opts.atol);

        if opts.fixed_step.is_none() && (err_norm > 1.0 || !err_norm.is_finite()) {
            let factor = if err_norm.is_finite() {
                MIN_FACTOR.max(SAFETY * err_norm.powf(-1.0 / (order as f64 + 1.0)))
            } else {
                MIN_FACTOR
            };
            change_step(&mut d, order, factor);
            h *= factor;
            n_equal = 0;
            stats.rejected += 1;
            continue;
        }

        // accept: update the differences
        stats.steps += 1;
        stats.steps_per_order[order] += 1;
        n_equal += 1;
        for e in 0..n {
            d[order + 2][e] = dcorr[e] - d[order + 1][e];
            d[order + 1][e] = dcorr[e];
        }
        for i in (0..=order).rev() {
            let (lo, hi) = d.split_at_mut(i + 1);
            for (a, b) in lo[i].iter_mut().zip(&hi[0]) {
                *a += *b;
            }
        }
        // output
        if t_eval.is_empty() {
            out_t.push(t_new);
            out_y.push(d[0].clone());
        } else {
            while next_eval < t_eval.len() && t_eval[next_eval] <= t_new + 1e-13 * span {
                let x = (t_eval[next_eval] - t_new) / h;
                let mut v = vec![T::ZERO; n];
                history_value(&d, order, x.min(0.0), &mut v);
                out_t.push(t_eval[next_eval]);
                out_y.push(v);
                next_eval += 1;
            }
        }
        t = t_new;
        if t >= t_end - 1e-13 * span {
            break;
        }

        if opts.fixed_step.is_some() {
            // only the order ramps up (no step changes)
            if n_equal > order && order < opts.max_order {
                order += 1;
                n_equal = 0;
            }
            continue;
        }
        if n_equal < order + 1 {
            continue;
        }
        // order / step selection
        let scale_y = &d[0];
        let en = |coef: f64, row: &[T]| -> f64 {
            let e: Vec<T> = row.iter().map(|v| v.scale(coef)).collect();
            scaled_rms(&e, scale_y, opts.rtol, opts.atol)
        };
        let err_m = if order > 1 { en(1.0 / order as f64, &d[order]) } else { f64::INFINITY };
        let err_p = if order < opts.max_order { en(1.0 / (order as f64 + 2.0), &d[order + 2]) } else { f64::INFINITY };
        let norms = [err_m, err_norm, err_p];
        let factors: Vec<f64> = norms
            .iter()
            .enumerate()
            .map(|(i, &e)| if e == 0.0 { f64::INFINITY } else { e.powf(-1.0 / (order as f64 + i as f64)) })
            .collect();
        let best = (0..3).max_by(|&a, &b| factors[a].total_cmp(&factors[b])).unwrap();
        let new_order = (order as isize + best as isize - 1) as usize;
        let factor = MAX_FACTOR.min(SAFETY * factors[best]);
        order = new_order.clamp(1, opts.max_order);
        if factor.is_finite() && !(0.9..1.2).contains(&factor) {
            let factor = factor.min(opts.max_step / h);
            change_step(&mut d, order, factor);
            h *= factor;
        }
        n_equal = 0;
    }

    stats.absorb(itm.fac.stats());
    stats.absorb(mfac.stats());
    stats.seconds = started.elapsed().as_secs_f64();
    Ok(Solution { t: out_t, y: out_y, stats })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_change_with_factor_one_is_the_identity() {
        let mut d: Vec<Vec<f64>> = (0..8).map(|i| vec![i as f64 + 1.0, (i * i) as f64 - 2.0]).collect();
        let before = d.clone();
        for order in 1..=5 {
            change_step(&mut d, order, 1.0);
            for i in 0..=order {
                for e in 0..2 {
                    assert!((d[i][e] - before[i][e]).abs() < 1e-12, "order {order} row {i}");
                }
            }
        }
    }

    #[test]
    fn step_change_preserves_polynomials_up_to_the_order() {
        // History of p(t) = 1 + 2t − t² + 0.5 t³ sampled with step h; after rescaling to step f·h the
        // differences must be those of the same polynomial on the new grid (anchored at the same t_n).
        let p = |t: f64| 1.0 + 2.0 * t - t * t + 0.5 * t * t * t;
        let (h, tn, order) = (0.1, 1.0, 3usize);
        let backward_differences = |h: f64| -> Vec<Vec<f64>> {
            let mut vals: Vec<f64> = (0..=order).map(|j| p(tn - j as f64 * h)).collect(); // y_n, y_{n-1}, ...
            let mut out = vec![vec![0.0]; 8];
            for k in 0..=order {
                out[k][0] = vals[0];
                for j in 0..order - k {
                    vals[j] -= vals[j + 1];
                }
            }
            out
        };
        for factor in [0.5, 1.7, 2.0] {
            let mut d = backward_differences(h);
            change_step(&mut d, order, factor);
            let want = backward_differences(h * factor);
            for k in 0..=order {
                assert!(
                    (d[k][0] - want[k][0]).abs() < 1e-10,
                    "factor {factor}, difference {k}: {} vs {}",
                    d[k][0],
                    want[k][0]
                );
            }
        }
    }
}
