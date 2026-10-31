//! Newmark-β and generalised-α (Chung–Hulbert 1993) for `M u'' + C u' + K u = f(t)`, fixed step.
//!
//! With `α_m = α_f = 0` this is Newmark-β. The default `β = 1/4, γ = 1/2` (average acceleration,
//! trapezoidal rule) is unconditionally stable, second-order accurate and has *no* numerical
//! dissipation: for `C = 0, f = 0` the discrete energy `½ vᵀMv + ½ uᵀKu` is conserved to rounding.
//! Generalised-α with spectral radius `ρ∞ ∈ [0, 1]` at infinite frequency damps the spurious
//! high-frequency modes of a spatial discretisation (`ρ∞ = 1` is Newmark average acceleration) and
//! keeps second-order accuracy:
//!
//! ```text
//! α_m = (2ρ∞ − 1)/(ρ∞ + 1),  α_f = ρ∞/(ρ∞ + 1),  γ = ½ − α_m + α_f,  β = ¼ (1 − α_m + α_f)²
//! ```
//!
//! One step solves `S a_{n+1} = rhs` with the constant matrix
//! `S = (1 − α_m) M + (1 − α_f) γ h C + (1 − α_f) β h² K`, factorised once for the whole run.

use super::{TsError, TsStats};
use crate::solvers::csr::Csr;
use crate::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::solvers::faer_glue::FaerScalar;
use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NewmarkParams {
    pub alpha_m: f64,
    pub alpha_f: f64,
    pub beta: f64,
    pub gamma: f64,
}

impl NewmarkParams {
    /// Newmark-β (`α_m = α_f = 0`).
    pub fn newmark(beta: f64, gamma: f64) -> Self {
        Self { alpha_m: 0.0, alpha_f: 0.0, beta, gamma }
    }
    /// `β = 1/4, γ = 1/2`: second order, energy conserving for the undamped homogeneous problem.
    pub fn average_acceleration() -> Self {
        Self::newmark(0.25, 0.5)
    }
    /// Generalised-α with high-frequency spectral radius `ρ∞ ∈ [0, 1]`.
    pub fn generalized_alpha(rho_inf: f64) -> Self {
        let alpha_m = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
        let alpha_f = rho_inf / (rho_inf + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f).powi(2);
        Self { alpha_m, alpha_f, beta, gamma }
    }
}

/// States at the output times.
#[derive(Debug, Clone)]
pub struct SecondOrderSolution<T> {
    pub t: Vec<f64>,
    pub u: Vec<Vec<T>>,
    pub v: Vec<Vec<T>>,
    pub stats: TsStats,
}

/// Takes `steps` steps of size `(t_end − t0)/steps`, storing the state every `output_every` steps
/// (and at `t0` and `t_end`).
#[allow(clippy::too_many_arguments)]
pub fn newmark_solve<T: FaerScalar>(
    m: &Csr<T>,
    c: Option<&Csr<T>>,
    k: &Csr<T>,
    force: &mut dyn FnMut(f64, &mut [T]),
    t0: f64,
    t_end: f64,
    steps: usize,
    u0: &[T],
    v0: &[T],
    params: NewmarkParams,
    output_every: usize,
    definiteness: Definiteness,
) -> Result<SecondOrderSolution<T>, TsError> {
    let n = u0.len();
    if v0.len() != n || m.nrows() != n || k.nrows() != n || c.is_some_and(|c| c.nrows() != n) {
        return Err(TsError::InvalidArgument("dimension mismatch between M, C, K, u0 and v0"));
    }
    if steps == 0 || !(t_end > t0) || output_every == 0 {
        return Err(TsError::InvalidArgument("need steps > 0, t_end > t0 and output_every > 0"));
    }
    let started = Instant::now();
    let h = (t_end - t0) / steps as f64;
    let NewmarkParams { alpha_m, alpha_f, beta, gamma } = params;
    // S = (1−αm) M + (1−αf) γ h C + (1−αf) β h² K
    let mut s_mat = Csr::axpby(T::from_f64(1.0 - alpha_m), m, T::from_f64((1.0 - alpha_f) * beta * h * h), k);
    if let Some(c) = c {
        s_mat = Csr::axpby(T::ONE, &s_mat, T::from_f64((1.0 - alpha_f) * gamma * h), c);
    }
    let mut sfac = ShiftedFactorization::new(&s_mat, None, definiteness)?;
    sfac.set_shift(T::ZERO)?;
    let mut mfac = ShiftedFactorization::new(m, None, definiteness)?;
    mfac.set_shift(T::ZERO)?;

    let mut stats = TsStats::default();
    let (mut u, mut v) = (u0.to_vec(), v0.to_vec());
    let mut f = vec![T::ZERO; n];
    let (mut tmp, mut tmp2) = (vec![T::ZERO; n], vec![T::ZERO; n]);
    // a0 = M⁻¹ (f(t0) − C v0 − K u0)
    force(t0, &mut f);
    stats.forcing_evals += 1;
    k.spmv(&u, &mut tmp);
    let mut a: Vec<T> = (0..n).map(|i| f[i] - tmp[i]).collect();
    if let Some(c) = c {
        c.spmv(&v, &mut tmp);
        for i in 0..n {
            a[i] -= tmp[i];
        }
    }
    mfac.solve_block_in_place(&mut a, 1);

    let mut out =
        SecondOrderSolution { t: vec![t0], u: vec![u.clone()], v: vec![v.clone()], stats: TsStats::default() };
    let (mut u_pred, mut v_pred, mut rhs) = (vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n]);
    for step in 0..steps {
        let tn = t0 + step as f64 * h;
        for i in 0..n {
            u_pred[i] = u[i] + v[i].scale(h) + a[i].scale(h * h * (0.5 - beta));
            v_pred[i] = v[i] + a[i].scale(h * (1.0 - gamma));
        }
        force(tn + (1.0 - alpha_f) * h, &mut f);
        stats.forcing_evals += 1;
        // rhs = f − αm M a − C[(1−αf) ṽ + αf v] − K[(1−αf) ũ + αf u]
        m.spmv(&a, &mut tmp);
        for i in 0..n {
            rhs[i] = f[i] - tmp[i].scale(alpha_m);
        }
        for i in 0..n {
            tmp2[i] = u_pred[i].scale(1.0 - alpha_f) + u[i].scale(alpha_f);
        }
        k.spmv(&tmp2, &mut tmp);
        for i in 0..n {
            rhs[i] -= tmp[i];
        }
        if let Some(c) = c {
            for i in 0..n {
                tmp2[i] = v_pred[i].scale(1.0 - alpha_f) + v[i].scale(alpha_f);
            }
            c.spmv(&tmp2, &mut tmp);
            for i in 0..n {
                rhs[i] -= tmp[i];
            }
        }
        sfac.solve_block_in_place(&mut rhs, 1);
        for i in 0..n {
            a[i] = rhs[i];
            u[i] = u_pred[i] + a[i].scale(beta * h * h);
            v[i] = v_pred[i] + a[i].scale(gamma * h);
        }
        stats.steps += 1;
        if (step + 1) % output_every == 0 || step + 1 == steps {
            out.t.push(t0 + (step + 1) as f64 * h);
            out.u.push(u.clone());
            out.v.push(v.clone());
        }
    }
    stats.absorb(sfac.stats());
    stats.absorb(mfac.stats());
    stats.seconds = started.elapsed().as_secs_f64();
    out.stats = stats;
    Ok(out)
}
