//! The θ-method for `M y' + K y = f(t)`, fixed step: `θ = 1` backward Euler (L-stable, first order),
//! `θ = ½` Crank–Nicolson, which for a linear problem is the **implicit midpoint rule** (A-stable,
//! second order, no damping of stiff modes: use BDF or generalised-α when those matter).
//!
//! `(M + θ h K) y_{n+1} = (M − (1 − θ) h K) y_n + h [θ f(t_{n+1}) + (1 − θ) f(t_n)]`; one
//! factorisation for the whole run.

use super::{Solution, TsError, TsStats};
use crate::linalg::csr::Csr;
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::linalg::faer_glue::FaerScalar;
use std::time::Instant;

#[allow(clippy::too_many_arguments)]
pub fn theta_solve<T: FaerScalar>(
    m: &Csr<T>,
    k: &Csr<T>,
    force: &mut dyn FnMut(f64, &mut [T]),
    t0: f64,
    t_end: f64,
    steps: usize,
    y0: &[T],
    theta: f64,
    output_every: usize,
    definiteness: Definiteness,
) -> Result<Solution<T>, TsError> {
    let n = y0.len();
    if m.nrows() != n || k.nrows() != n {
        return Err(TsError::InvalidArgument("dimension mismatch between M, K and y0"));
    }
    if steps == 0 || !(t_end > t0) || output_every == 0 || !(0.0..=1.0).contains(&theta) || theta == 0.0 {
        return Err(TsError::InvalidArgument("need steps > 0, t_end > t0, output_every > 0 and 0 < θ ≤ 1"));
    }
    let started = Instant::now();
    let h = (t_end - t0) / steps as f64;
    // M + θ h K = K − σ M with σ = −1/(θ h), scaled: use the explicit matrix instead.
    let lhs = Csr::axpby(T::ONE, m, T::from_f64(theta * h), k);
    let mut fac = ShiftedFactorization::new(&lhs, None, definiteness)?;
    fac.set_shift(T::ZERO)?;
    let mut stats = TsStats::default();
    let mut y = y0.to_vec();
    let (mut f0, mut f1) = (vec![T::ZERO; n], vec![T::ZERO; n]);
    let (mut tmp, mut rhs) = (vec![T::ZERO; n], vec![T::ZERO; n]);
    force(t0, &mut f0);
    stats.forcing_evals += 1;
    let mut out = Solution { t: vec![t0], y: vec![y.clone()], stats: TsStats::default() };
    for step in 0..steps {
        let t1 = t0 + (step + 1) as f64 * h;
        force(t1, &mut f1);
        stats.forcing_evals += 1;
        m.spmv(&y, &mut rhs);
        k.spmv(&y, &mut tmp);
        for i in 0..n {
            rhs[i] += -tmp[i].scale((1.0 - theta) * h) + (f1[i].scale(theta) + f0[i].scale(1.0 - theta)).scale(h);
        }
        fac.solve_block_in_place(&mut rhs, 1);
        y.copy_from_slice(&rhs);
        std::mem::swap(&mut f0, &mut f1);
        stats.steps += 1;
        if (step + 1) % output_every == 0 || step + 1 == steps {
            out.t.push(t1);
            out.y.push(y.clone());
        }
    }
    stats.absorb(fac.stats());
    stats.seconds = started.elapsed().as_secs_f64();
    out.stats = stats;
    Ok(out)
}
