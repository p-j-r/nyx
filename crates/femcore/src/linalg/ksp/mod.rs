//! Krylov linear solvers (the role of PETSc `KSP`) for real and complex systems.
//!
//! * [`Cg`]: preconditioned conjugate gradients; [`CgKind::Hermitian`] for Hermitian positive
//!   definite systems, [`CgKind::Symmetric`] = **COCG** (conjugate gradients with the bilinear
//!   form `xᵀy`, no conjugation) for complex *symmetric* systems such as PML problems. PETSc's
//!   `KSPCG` with `-ksp_cg_type symmetric` is the same algorithm.
//! * [`Minres`]: Hermitian *indefinite* systems with an SPD preconditioner; [`Bicgstab`]: short-recurrence
//!   solver for general (also complex) systems, right-preconditioned.
//! * [`Cocr`]: conjugate residuals for complex *symmetric* systems (smoother than COCG); [`Idrs`]: IDR(s).
//! * [`Qmr`]: quasi-minimal residual (smooth BiCG, two products per iteration, or one for complex
//!   symmetric systems); [`GcroDr`]: GMRES with Krylov-subspace recycling for sequences of systems.
//! * [`Gmres`]: restarted GMRES(m) with left or right preconditioning and the flexible variant
//!   (FGMRES, needed when the preconditioner changes from step to step, e.g. inner iterations or
//!   AMG with different cycles).
//!
//! Defaults, norms and the convergence test follow PETSc 3.19 (`rtol 1e-5`, `atol 1e-50`,
//! `dtol 1e4`, `max_it 10000`; CG monitors the preconditioned residual norm; GMRES restarts
//! every 30 iterations with left preconditioning and classical Gram–Schmidt), so iteration
//! counts can be compared with PETSc one to one (`harness/`).

mod bicgstab;
mod cg;
mod cocr;
mod gcrodr;
mod gmres;
mod idrs;
mod minres;
mod qmr;

pub use bicgstab::{Bicgstab, BicgstabOptions};
pub use cg::{Cg, CgKind, CgNorm, CgOptions};
pub use cocr::{Cocr, CocrOptions};
pub use gcrodr::{GcroDr, GcroDrOptions, GcroDrStats};
pub use gmres::{Gmres, GmresOptions, Orthogonalization, PcSide};
pub use idrs::{Idrs, IdrsOptions};
pub use minres::{Minres, MinresOptions};
pub use qmr::{Qmr, QmrOptions};

use crate::linalg::scalar::Scalar;

/// Why a solve stopped (PETSc `KSPConvergedReason`, the cases that can occur here).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvergedReason {
    /// `‖r‖ ≤ rtol · ‖r₀‖`.
    Rtol,
    /// `‖r‖ < atol`.
    Atol,
    DivergedIts,
    /// `‖r‖ ≥ dtol · ‖r₀‖`.
    DivergedDtol,
    /// Breakdown (GMRES, COCG `pᵀAp = 0`) or non-finite scalar.
    DivergedBreakdown,
    /// Hermitian CG saw `pᴴAp ≤ 0`: `A` is not positive definite.
    DivergedIndefiniteMat,
    /// Hermitian CG saw `rᴴM⁻¹r` change sign: `M` is not positive definite.
    DivergedIndefinitePc,
    DivergedNanOrInf,
}

impl ConvergedReason {
    pub fn is_converged(self) -> bool {
        matches!(self, ConvergedReason::Rtol | ConvergedReason::Atol)
    }
}

/// Result of a linear solve.
#[derive(Debug, Clone)]
pub struct SolveInfo {
    pub iterations: usize,
    pub reason: ConvergedReason,
    /// The monitored residual norm at exit (see each solver's norm/side option).
    pub residual_norm: f64,
    /// One entry per iteration (index 0 is the initial norm) if `record_history` is set: the
    /// convergence-history log used to compare with PETSc's `-ksp_monitor`.
    pub history: Vec<f64>,
    pub matvecs: usize,
    pub pc_applies: usize,
}

/// Convergence tolerances (PETSc `KSPSetTolerances`).
#[derive(Debug, Clone, Copy)]
pub struct Tolerances {
    pub rtol: f64,
    pub atol: f64,
    pub dtol: f64,
    pub max_it: usize,
}

impl Default for Tolerances {
    fn default() -> Self {
        Self { rtol: 1e-5, atol: 1e-50, dtol: 1e4, max_it: 10_000 }
    }
}

/// PETSc `KSPConvergedDefault`.
#[derive(Debug, Clone)]
pub(crate) struct ConvTest {
    tol: Tolerances,
    rnorm0: f64,
    ttol: f64,
}

impl ConvTest {
    pub(crate) fn new(tol: Tolerances) -> Self {
        Self { tol, rnorm0: 0.0, ttol: 0.0 }
    }
    pub(crate) fn set_reference(&mut self, rnorm0: f64) {
        self.rnorm0 = rnorm0;
        self.ttol = (self.tol.rtol * rnorm0).max(self.tol.atol);
    }
    pub(crate) fn rnorm0(&self) -> f64 {
        self.rnorm0
    }
    pub(crate) fn check(&self, rnorm: f64) -> Option<ConvergedReason> {
        if !rnorm.is_finite() {
            Some(ConvergedReason::DivergedNanOrInf)
        } else if rnorm <= self.ttol {
            Some(if rnorm < self.tol.atol { ConvergedReason::Atol } else { ConvergedReason::Rtol })
        } else if rnorm >= self.tol.dtol * self.rnorm0 {
            Some(ConvergedReason::DivergedDtol)
        } else {
            None
        }
    }
}

pub(crate) fn ensure_len<T: Scalar>(v: &mut Vec<T>, n: usize) {
    if v.len() != n {
        v.clear();
        v.resize(n, T::ZERO);
    }
}
