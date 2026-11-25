//! Rosenbrock (ROS2, Verwer–Spee–Blom–Hundsdorfer 1999) for the stiff nonlinear autonomous problem
//! `M y' = F(y)` with a user-supplied Jacobian `J = ∂F/∂y` (sparse), e.g. heat conduction with
//! temperature-dependent conductivity. Time-dependent problems: append `t` as an unknown
//! (`t' = 1`).
//!
//! `γ = 1 + 1/√2`, second order, L-stable, independent of the Jacobian's accuracy for stability:
//!
//! ```text
//! (M − γ h J) k₁ = F(y_n)
//! (M − γ h J) k₂ = F(y_n + h k₁) − 2 M k₁
//! y_{n+1} = y_n + h (3/2 k₁ + 1/2 k₂)
//! ```
//!
//! Error control by **step doubling**: one step of `h` against two of `h/2` (Jacobians at `y_n` and at
//! the midpoint); the doubled result is kept and `‖y_{h/2,2} − y_h‖ / 3` is the local error estimate,
//! `h_new = h · clamp(0.9 · err^{-1/3}, 0.2, 5)`. Three numeric factorisations per step, all on one
//! symbolic analysis per Jacobian object (values are replaced with
//! [`ShiftedFactorization::replace_matrices`]); an embedded estimator would cut that, see
//! `docs/STATUS.md`.

use super::{Solution, TsError, TsStats, scaled_rms};
use crate::linalg::csr::Csr;
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::linalg::faer_glue::FaerScalar;
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct RosenbrockOptions {
    pub rtol: f64,
    pub atol: f64,
    pub first_step: f64,
    pub max_step: f64,
    /// Constant step, no error control (for convergence studies).
    pub fixed_step: bool,
    pub max_steps: usize,
    pub definiteness: Definiteness,
}

impl RosenbrockOptions {
    pub fn new(rtol: f64, atol: f64, first_step: f64) -> Self {
        Self {
            rtol,
            atol,
            first_step,
            max_step: f64::INFINITY,
            fixed_step: false,
            max_steps: 1_000_000,
            definiteness: Definiteness::General,
        }
    }
}

const GAMMA: f64 = 1.0 + std::f64::consts::FRAC_1_SQRT_2;

struct Stage<'a, T: FaerScalar> {
    m: &'a Csr<T>,
    f: &'a mut dyn FnMut(&[T], &mut [T]),
    stats: &'a mut TsStats,
    factor_secs: f64,
    numeric: usize,
    symbolic: usize,
}

impl<T: FaerScalar> Stage<'_, T> {
    /// One ROS2 step of size `h` from `y` with the Jacobian `j`, through a factorisation of
    /// `M − γ h J` (= `A − σ B` with `A = M`, `B = J`, `σ = γ h`).
    fn step(&mut self, fac: &mut ShiftedFactorization<T>, h: f64, y: &[T], out: &mut [T]) -> Result<(), TsError> {
        let n = y.len();
        fac.set_shift(T::from_f64(GAMMA * h))?;
        let mut k1 = vec![T::ZERO; n];
        (self.f)(y, &mut k1);
        self.stats.forcing_evals += 1;
        fac.solve_block_in_place(&mut k1, 1)?;
        let y2: Vec<T> = y.iter().zip(&k1).map(|(&a, &b)| a + b.scale(h)).collect();
        let mut k2 = vec![T::ZERO; n];
        (self.f)(&y2, &mut k2);
        self.stats.forcing_evals += 1;
        let mut mk1 = vec![T::ZERO; n];
        self.m.spmv(&k1, &mut mk1);
        for i in 0..n {
            k2[i] -= mk1[i].scale(2.0);
        }
        fac.solve_block_in_place(&mut k2, 1)?;
        for i in 0..n {
            out[i] = y[i] + (k1[i].scale(1.5) + k2[i].scale(0.5)).scale(h);
        }
        Ok(())
    }
    fn account(&mut self, fac: &ShiftedFactorization<T>) {
        let s = fac.stats();
        self.numeric += s.numeric_count;
        self.symbolic += s.symbolic_count;
        self.factor_secs += s.numeric_secs_total;
    }
}

/// Integrates `M y' = F(y)` from `t0` to `t_end`; output at `t_eval` by linear interpolation of the
/// accepted steps (empty = every accepted step). `jac(y)` returns `∂F/∂y` at `y`; all Jacobians must
/// share one sparsity pattern with `M`.
#[allow(clippy::too_many_arguments)]
pub fn rosenbrock_solve<T: FaerScalar>(
    m: &Csr<T>,
    f: &mut dyn FnMut(&[T], &mut [T]),
    jac: &mut dyn FnMut(&[T]) -> Csr<T>,
    t0: f64,
    t_end: f64,
    y0: &[T],
    t_eval: &[f64],
    opts: &RosenbrockOptions,
) -> Result<Solution<T>, TsError> {
    let n = y0.len();
    if m.nrows() != n || !(t_end > t0) || !(opts.first_step > 0.0) {
        return Err(TsError::InvalidArgument("need M n×n, t_end > t0 and first_step > 0"));
    }
    let started = Instant::now();
    let mut stats = TsStats::default();
    let mut st = Stage { m, f, stats: &mut stats, factor_secs: 0.0, numeric: 0, symbolic: 0 };
    let (mut t, mut h) = (t0, opts.first_step.min(opts.max_step));
    let mut y = y0.to_vec();
    let span = t_end - t0;
    let mut out_t = vec![t0];
    let mut out_y = vec![y.clone()];
    let mut taken = 0usize;
    // two persistent factorisations: values are replaced every step, the symbolic analysis is kept
    let (mut fac0, mut fac1): (Option<ShiftedFactorization<T>>, Option<ShiftedFactorization<T>>) = (None, None);
    let (mut y_full, mut y_half, mut y_two) = (vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n]);
    while t < t_end - 1e-13 * span {
        taken += 1;
        if taken > opts.max_steps {
            return Err(TsError::TooManySteps(t));
        }
        h = h.min(t_end - t).min(opts.max_step);
        if h < 1e-14 * span.max(1.0) {
            return Err(TsError::StepTooSmall(t));
        }
        let j0 = jac(&y);
        let fac0 = match fac0.as_mut() {
            Some(f) => {
                f.replace_matrices(m, Some(&j0))?;
                f
            }
            None => fac0.insert(ShiftedFactorization::new(m, Some(&j0), opts.definiteness)?),
        };
        if opts.fixed_step {
            st.step(fac0, h, &y, &mut y_full)?;
            y.copy_from_slice(&y_full);
            t += h;
            st.stats.steps += 1;
            out_t.push(t);
            out_y.push(y.clone());
            continue;
        }
        st.step(fac0, h, &y, &mut y_full)?;
        st.step(fac0, 0.5 * h, &y, &mut y_half)?;
        let j1 = jac(&y_half);
        let fac1 = match fac1.as_mut() {
            Some(f) => {
                f.replace_matrices(m, Some(&j1))?;
                f
            }
            None => fac1.insert(ShiftedFactorization::new(m, Some(&j1), opts.definiteness)?),
        };
        st.step(fac1, 0.5 * h, &y_half, &mut y_two)?;
        let diff: Vec<T> = y_two.iter().zip(&y_full).map(|(&a, &b)| (a - b).scale(1.0 / 3.0)).collect();
        let err = scaled_rms(&diff, &y_two, opts.rtol, opts.atol);
        if err.is_finite() && err <= 1.0 {
            t += h;
            y.copy_from_slice(&y_two);
            st.stats.steps += 1;
            out_t.push(t);
            out_y.push(y.clone());
        } else {
            st.stats.rejected += 1;
        }
        let factor = if err.is_finite() { (0.9 * err.max(1e-10).powf(-1.0 / 3.0)).clamp(0.2, 5.0) } else { 0.2 };
        h *= factor;
    }
    for f in [&fac0, &fac1].into_iter().flatten() {
        st.account(f);
    }
    st.stats.factorizations = st.numeric;
    st.stats.symbolic_analyses = st.symbolic;
    st.stats.factor_seconds = st.factor_secs;
    stats.seconds = started.elapsed().as_secs_f64();
    if t_eval.is_empty() {
        return Ok(Solution { t: out_t, y: out_y, stats });
    }
    // linear interpolation between accepted steps
    let mut ys = Vec::with_capacity(t_eval.len());
    for &te in t_eval {
        let i = out_t.partition_point(|&x| x < te).clamp(1, out_t.len() - 1);
        let (ta, tb) = (out_t[i - 1], out_t[i]);
        let w = if tb > ta { ((te - ta) / (tb - ta)).clamp(0.0, 1.0) } else { 1.0 };
        ys.push(out_y[i - 1].iter().zip(&out_y[i]).map(|(&a, &b)| a.scale(1.0 - w) + b.scale(w)).collect());
    }
    Ok(Solution { t: t_eval.to_vec(), y: ys, stats })
}
