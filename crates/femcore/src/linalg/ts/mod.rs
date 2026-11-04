//! Time integration of linear semi-discrete systems (port item T1).
//!
//! * [`bdf`]: **variable-step, variable-order BDF** (orders 1–5) with error control, for the
//!   first-order system `M y' + K y = f(t)` (heat, diffusion, damped first-order wave forms).
//! * [`newmark`]: **Newmark-β and generalised-α** (Chung–Hulbert) for the second-order system
//!   `M u'' + C u' + K u = f(t)` (elastic waves); fixed step, one factorisation for the whole run.
//! * [`theta`]: the **θ-method** (backward Euler, Crank–Nicolson = implicit midpoint for linear
//!   problems) for the first-order system, fixed step.
//! * [`rosenbrock`]: **Rosenbrock (ROS2)** for the nonlinear `M y' = F(t, y)` with a user Jacobian
//!   (stiff nonlinear heat), step-doubling error control.
//!
//! All use [`crate::linalg::direct::ShiftedFactorization`] so the symbolic analysis is shared by every step size, and a
//! numeric factorisation is reused for as long as the step size and order do not change. Tolerances
//! are explicit arguments (no hidden defaults): the port list warns that other codes' default
//! tolerances were inaccurate.
//!
//! The solvers are generic over real and complex scalars; time is always `f64`.

#![allow(clippy::needless_range_loop, clippy::neg_cmp_op_on_partial_ord)] // index loops over several vectors; `!(x > 0.0)` is deliberate (rejects NaN)

pub mod bdf;
pub mod newmark;
pub mod rosenbrock;
pub mod theta;

use crate::linalg::direct::{FactorStats, SolverError};
use crate::linalg::scalar::Scalar;

#[derive(Debug, Clone, PartialEq)]
pub enum TsError {
    Solver(SolverError),
    InvalidArgument(&'static str),
    /// The step size fell below `min_step` at this time.
    StepTooSmall(f64),
    /// `max_steps` accepted + rejected steps were taken before reaching the end time.
    TooManySteps(f64),
}

impl std::fmt::Display for TsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TsError::Solver(e) => write!(f, "linear solve failed: {e}"),
            TsError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
            TsError::StepTooSmall(t) => write!(f, "step size underflow at t = {t}"),
            TsError::TooManySteps(t) => write!(f, "step limit reached at t = {t}"),
        }
    }
}

impl std::error::Error for TsError {}

impl From<SolverError> for TsError {
    fn from(e: SolverError) -> Self {
        TsError::Solver(e)
    }
}

/// Counters of one integration.
#[derive(Debug, Clone, Default)]
pub struct TsStats {
    pub steps: usize,
    pub rejected: usize,
    /// Numeric factorisations and symbolic analyses performed.
    pub factorizations: usize,
    pub symbolic_analyses: usize,
    pub solves: usize,
    pub forcing_evals: usize,
    /// Accepted steps per BDF order (index = order; unused for the fixed-step methods).
    pub steps_per_order: [usize; 6],
    pub seconds: f64,
    pub factor_seconds: f64,
    pub solve_seconds: f64,
}

impl TsStats {
    /// Adds the counters of one factorisation object (call once per object at the end of a run).
    pub(crate) fn absorb(&mut self, f: FactorStats) {
        self.factorizations += f.numeric_count;
        self.symbolic_analyses += f.symbolic_count;
        self.factor_seconds += f.numeric_secs_total;
        self.solve_seconds += f.solve_secs_total;
        self.solves += f.solves;
    }
}

/// Output of an integration: the state at the requested times.
#[derive(Debug, Clone)]
pub struct Solution<T> {
    /// Output times (the requested `t_eval`, or the accepted step times when none were requested).
    pub t: Vec<f64>,
    pub y: Vec<Vec<T>>,
    pub stats: TsStats,
}

/// Root-mean-square norm of `e_i / (atol + rtol·|y_i|)`.
pub(crate) fn scaled_rms<T: Scalar>(e: &[T], y: &[T], rtol: f64, atol: f64) -> f64 {
    let s: f64 = e.iter().zip(y).map(|(&ei, &yi)| ei.modulus_sq() / (atol + rtol * yi.modulus()).powi(2)).sum();
    (s / e.len().max(1) as f64).sqrt()
}
