//! Eigenvalue solvers (the role of SLEPc `EPS`).
//!
//! [`krylov_schur`] is a thick-restart Krylov–Schur method with locking, written for real and
//! complex scalars, with an optional `B`-weighted inner product and an optional block size
//! (block variants find repeated eigenvalues, which single-vector Krylov methods do not do
//! reliably).
//!
//! Two projected-problem kinds are supported:
//! * [`Kind::Hermitian`]: the operator is self-adjoint in the inner product used (real
//!   symmetric, complex Hermitian, or `B`-self-adjoint). The projected problem is solved by a
//!   Hermitian eigendecomposition and the Schur vectors are the Ritz vectors.
//! * [`Kind::General`]: any complex operator (complex symmetric PML problems, for instance),
//!   using a complex Schur form with eigenvalue reordering. Real non-symmetric operators are
//!   not supported (their Schur form has 2×2 blocks); promote them to complex.

mod krylov_schur;
mod orth;

pub use krylov_schur::krylov_schur;

use crate::linalg::basis::Basis;
use crate::linalg::dense::DenseError;
use crate::linalg::scalar::{C64, Scalar};

/// Which eigenvalues of the operator handed to the solver to compute (SLEPc `EPSWhich`).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Which {
    /// Largest `|λ|` (SLEPc default); with shift-and-invert this selects eigenvalues closest to
    /// the shift.
    #[default]
    LargestMagnitude,
    SmallestMagnitude,
    LargestReal,
    SmallestReal,
    LargestImag,
    SmallestImag,
    /// Closest to `τ` in the complex plane.
    TargetMagnitude(C64),
    /// Real part closest to `τ`.
    TargetReal(f64),
}

impl Which {
    /// Sort key: smaller means "more wanted".
    pub fn key(self, lambda: C64) -> f64 {
        match self {
            Which::LargestMagnitude => -lambda.norm(),
            Which::SmallestMagnitude => lambda.norm(),
            Which::LargestReal => -lambda.re,
            Which::SmallestReal => lambda.re,
            Which::LargestImag => -lambda.im.abs(),
            Which::SmallestImag => lambda.im.abs(),
            Which::TargetMagnitude(t) => (lambda - t).norm(),
            Which::TargetReal(t) => (lambda.re - t).abs(),
        }
    }
}

/// How the projected problem is solved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Hermitian,
    General,
}

#[derive(Debug, Clone)]
pub struct EpsOptions {
    /// Number of eigenpairs wanted.
    pub nev: usize,
    /// Maximum dimension of the Krylov basis; default `max(2·nev, nev + 15)` (SLEPc), raised to
    /// leave room for the block size, and capped by the problem size.
    pub ncv: Option<usize>,
    /// Relative residual tolerance `‖Ax − λx‖ ≤ tol·|λ|` for the estimated residual (SLEPc
    /// default `1e-8`).
    pub tol: f64,
    /// Maximum number of restarts; default `max(100, 2n/ncv)`.
    pub max_it: Option<usize>,
    pub which: Which,
    /// Fraction of the unconverged Ritz vectors kept at each restart (SLEPc default `0.5`).
    pub keep: f64,
    /// Vectors per block; `1` is the standard single-vector method.
    pub block_size: usize,
    pub kind: Kind,
    /// Measure convergence on the *original* problem of a shift-and-invert transformation
    /// `λ = σ + 1/θ`: `Some(σ)` tests `‖r‖/(|θ|·|λ|)`, the relative residual `‖Ax−λBx‖/|λ|` of
    /// the original pencil, instead of SLEPc's `‖r‖/|θ|` (`None`, the default here so that
    /// iteration counts can be compared with SLEPc one to one).
    ///
    /// The two differ a lot when `|λ| ≪ |θ|⁻¹`, e.g. a shift far below small eigenvalues:
    /// `‖r‖/|θ| < tol` then only guarantees `‖Ax−λBx‖/|λ| < tol/|λ|`.
    /// [`crate::linalg::problem::KrylovSchurShiftInvert`] sets this automatically.
    pub original_shift: Option<C64>,
    /// Replace the start block `X₀` by `Op·X₀` before orthonormalising, forcing it into the range
    /// of the operator. SLEPc does this for generalised Hermitian problems (`EPSGetStartVector`):
    /// it removes components in the null space of a singular `B` (infinite eigenvalues), which
    /// otherwise pollute the Krylov space. `None` = on when a weight matrix is given and the kind
    /// is Hermitian, off otherwise (SLEPc's rule); costs one extra operator application per
    /// start vector.
    pub start_in_range: Option<bool>,
    /// Seed of the pseudo-random vectors (extra start-block columns, breakdown replacements).
    pub seed: u64,
    /// Record per-restart diagnostics in [`EigenResult::history`].
    pub record_history: bool,
}

impl EpsOptions {
    fn base(nev: usize, kind: Kind) -> Self {
        Self {
            nev,
            ncv: None,
            tol: 1e-8,
            max_it: None,
            which: Which::default(),
            keep: 0.5,
            block_size: 1,
            kind,
            original_shift: None,
            start_in_range: None,
            seed: 0x5EED_0001,
            record_history: false,
        }
    }
    /// Options for a Hermitian (self-adjoint) operator.
    pub fn hermitian(nev: usize) -> Self {
        Self::base(nev, Kind::Hermitian)
    }
    /// Options for a general complex operator.
    pub fn general(nev: usize) -> Self {
        Self::base(nev, Kind::General)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EpsReason {
    /// At least `nev` eigenpairs converged.
    Converged,
    /// Restart limit reached.
    DivergedIts,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EpsError {
    InvalidArgument(&'static str),
    Unsupported(&'static str),
    Dense(DenseError),
}

impl std::fmt::Display for EpsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EpsError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
            EpsError::Unsupported(m) => write!(f, "unsupported: {m}"),
            EpsError::Dense(e) => write!(f, "projected problem failed: {e}"),
        }
    }
}

impl std::error::Error for EpsError {}

impl From<DenseError> for EpsError {
    fn from(e: DenseError) -> Self {
        EpsError::Dense(e)
    }
}

/// Seconds spent in each phase (the per-phase timers every solver here carries, so that a
/// slowdown can be located without a profiler).
#[derive(Debug, Clone, Copy, Default)]
pub struct Timings {
    /// Operator applications (for shift-and-invert: the factorisation solves).
    pub operator: f64,
    /// Products with the weight matrix `B` inside orthogonalisation.
    pub weight: f64,
    /// Orthogonalisation (dots and updates), excluding `B` products.
    pub orthogonalize: f64,
    /// Dense projected problems (Schur / Hermitian eigendecomposition, reordering).
    pub dense: f64,
    /// Basis updates at restarts (`V ← V·Q`).
    pub restart: f64,
    pub total: f64,
}

/// One entry per restart cycle when `record_history` is set (residual-history log, the
/// analogue of SLEPc's `-eps_monitor`).
#[derive(Debug, Clone)]
pub struct HistoryEntry {
    /// 1-based cycle number.
    pub cycle: usize,
    pub nconv: usize,
    /// Dimension of the projected problem at the end of the cycle.
    pub nv: usize,
    pub matvecs: usize,
    /// Relative error estimates of the leading Ritz pairs of the active window, in wanted order.
    pub errors: Vec<f64>,
    /// The corresponding Ritz values (of the operator given to the solver).
    pub ritz: Vec<C64>,
}

/// Converged eigenpairs, ordered from most to least wanted.
#[derive(Debug)]
pub struct EigenResult<T: Scalar> {
    pub reason: EpsReason,
    /// Eigenvalues of the operator given to the solver (complex; imaginary parts are zero for
    /// Hermitian problems).
    pub values: Vec<C64>,
    /// Relative error estimates `‖r‖/|λ|` for each returned pair.
    pub errors: Vec<f64>,
    /// Eigenvectors as the leading `values.len()` columns. Hermitian kind: orthonormal in the
    /// inner product used (`xᴴBy`, i.e. mass-normalised modes, or Euclidean without a weight);
    /// general kind: unit 2-norm.
    pub vectors: Basis<T>,
    /// Number of converged pairs (may exceed `nev`).
    pub nconv: usize,
    /// Krylov dimension actually used (after defaults and capping).
    pub ncv: usize,
    /// Number of restart cycles performed (SLEPc `EPSGetIterationNumber`).
    pub iterations: usize,
    /// Operator applications, counted per vector.
    pub matvecs: usize,
    /// Vectors multiplied by the weight matrix `B`.
    pub b_products: usize,
    /// Breakdowns (invariant subspaces found and replaced by fresh vectors).
    pub breakdowns: usize,
    /// Orthogonalisation calls and how many needed the DGKS refinement pass.
    pub orth_calls: usize,
    pub orth_refinements: usize,
    pub timings: Timings,
    pub history: Vec<HistoryEntry>,
}

impl<T: Scalar> EigenResult<T> {
    pub fn vector(&self, i: usize) -> &[T] {
        self.vectors.col(i)
    }
}
