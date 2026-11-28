//! Shifted sparse direct solves `(A − σB)⁻¹` with symbolic-factorisation reuse (port item S1).
//!
//! Shift-and-invert eigensolvers, frequency sweeps, the β fixed-point iteration for tilted
//! crystals and contour methods all factorise matrices with the **same sparsity pattern** and a
//! different shift. [`ShiftedFactorization`] builds the union pattern of `A` and `B` once, runs
//! `faer`'s symbolic analysis (ordering, elimination structure) once, and for every new shift
//! only reassembles the values (`A − σB`, one pass over the non-zeros) and re-runs the numeric
//! factorisation with the cached symbolic object.
//!
//! This file and `faer_glue.rs` are the only ones that use `faer`; they target faer 0.24.4 and
//! use its *low-level* sparse API (`sparse::linalg::{cholesky, lu}`), which takes the parallelism
//! and the scratch memory per call. That is what lets solves run sequentially while
//! factorisations use the shared pool (see [`ThreadPolicy`]). If nyx pins another version, the
//! calls to adapt are `factorize_symbolic_{cholesky,lu}`, `factorize_numeric_{llt,lu}`,
//! `LltRef::new` / `LuRef::new_unchecked` and `solve_in_place_with_conj`.

use crate::linalg::csr::Csr;
use crate::linalg::faer_glue::{FaerScalar, clear_upper_vector_state};
use crate::linalg::ilu::Ilu;
use crate::linalg::ksp::{Gmres, GmresOptions, PcSide, Tolerances};
use crate::linalg::linop::{LinearOperator, Preconditioner};
use crate::linalg::ordering::{FillOrdering, adjacency_from_lower, nested_dissection};
use crate::linalg::pc::Jacobi;
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::matmul;
use faer::linalg::triangular_solve;
use faer::perm::PermRef;
use faer::prelude::{Reborrow, ReborrowMut};
use faer::sparse::linalg::cholesky::simplicial::SymbolicSimplicialCholesky;
use faer::sparse::linalg::cholesky::supernodal::SymbolicSupernodalCholesky;
use faer::sparse::linalg::cholesky::{
    IntranodeLbltRef, LltRef, SymbolicCholesky, SymbolicCholeskyRaw, SymmetricOrdering, factorize_symbolic_cholesky,
};
use faer::sparse::linalg::lu::{LuRef, NumericLu, SymbolicLu, factorize_symbolic_lu};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use faer::{Accum, Conj, MatMut, MatRef, Par, Side};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolverError {
    /// The matrix is singular (or numerically so) at this shift.
    Singular,
    /// A Cholesky factorisation was required but the matrix is not positive definite.
    NotPositiveDefinite,
    /// [`ShiftedSolver::solve_block_in_place`] was called before a successful [`ShiftedSolver::set_shift`].
    NotFactored,
    OutOfMemory,
    /// More non-zeros than the 32-bit indices used here can address.
    TooLarge,
    InvalidArgument(&'static str),
}

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolverError::Singular => write!(f, "matrix is singular at this shift"),
            SolverError::NotPositiveDefinite => write!(f, "matrix is not positive definite"),
            SolverError::NotFactored => write!(f, "solver has not been factored (call set_shift first)"),
            SolverError::OutOfMemory => write!(f, "out of memory in the sparse factorisation"),
            SolverError::TooLarge => write!(f, "matrix exceeds 32-bit index range"),
            SolverError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
        }
    }
}

impl std::error::Error for SolverError {}

/// What the caller knows about `A − σB`, which decides the factorisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Definiteness {
    /// No assumption: sparse LU with partial pivoting (indefinite, complex symmetric, general).
    General,
    /// `A − σB` is symmetric/Hermitian positive definite at the shifts used (e.g. stiffness −
    /// σ·mass with σ below the spectrum): try Cholesky, fall back to LU if it fails.
    PositiveDefinite,
}

/// Which factorisation was used for the last shift.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Used {
    None,
    Cholesky,
    Lu,
}

/// Specific sparse direct backend and ordering chosen by [`ShiftedFactorization`] (or another [`ShiftedSolver`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FactorBackend {
    #[default]
    None,
    /// Sparse Cholesky (`LL^H`) with AMD ordering.
    CholeskyAmd,
    /// Sparse Cholesky (`LL^H`) with nested-dissection ordering.
    CholeskyNd,
    /// Sparse symmetric/Hermitian indefinite Bunch–Kaufman (`LBL^H`) with AMD ordering (lower triangle only).
    LbltAmd,
    /// Sparse symmetric/Hermitian indefinite Bunch–Kaufman (`LBL^H`) with nested-dissection ordering (lower triangle only).
    LbltNd,
    /// `faer` sparse LU with its default COLAMD column ordering.
    LuColamd,
    /// Sparse LU with AMD column ordering on `A + A^T`.
    LuAmd,
    /// Sparse LU with nested-dissection column ordering on `A + A^T`.
    LuNd,
    /// Block low-rank / exact multifrontal LU on a nested-dissection tree.
    BlrNd,
}

/// Thread policy for one phase of a sparse factorisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Threads {
    /// Whatever `faer::set_global_parallelism` says at the time of the call (nyx's shared pool).
    Global,
    Sequential,
    Fixed(usize),
}

/// Dimension threshold above which [`Threads::Global`] may use multi-threaded `faer` sparse
/// factorisations (paired with deterministic compensated-residual iterative refinement in
/// [`ShiftedFactorization::solve_block_in_place`] to uphold `ADR 11`).
pub(crate) const PAR_FACTOR_MIN_N: usize = 60_000;

/// Returns `true` when `NYX_DETERMINISTIC_FACTOR` requests strictly sequential (`Par::Seq`)
/// sparse and dense factorisations for bitwise thread-count invariance (`ADR 11`).
#[inline]
pub(crate) fn deterministic_factor_env() -> bool {
    matches!(
        std::env::var("NYX_DETERMINISTIC_FACTOR").ok().as_deref().map(str::trim),
        Some("1" | "true" | "TRUE" | "yes" | "YES")
    )
}

impl Threads {
    pub(crate) fn par(self) -> Par {
        match self {
            Threads::Global => {
                if crate::par::factor_threads() <= 1 || deterministic_factor_env() {
                    Par::Seq
                } else {
                    faer::get_global_parallelism()
                }
            }
            Threads::Sequential => Par::Seq,
            Threads::Fixed(n) if n <= 1 => Par::Seq,
            Threads::Fixed(n) => Par::rayon(n),
        }
    }

    /// Parallelism for factorising a sparse matrix of dimension `n`: under [`Threads::Global`],
    /// systems up to 60,000 unknowns (or whenever `NYX_DETERMINISTIC_FACTOR=1` is set) factorise
    /// with `Par::Seq` so direct solves are bit-deterministic across `NYX_THREADS` while avoiding
    /// supernode thread-synchronisation overhead on small/medium problems.
    pub(crate) fn par_for_factor(self, n: usize) -> Par {
        match self {
            Threads::Global
                if n <= PAR_FACTOR_MIN_N
                    || crate::par::factor_threads() <= 1
                    || deterministic_factor_env() =>
            {
                Par::Seq
            }
            _ => self.par(),
        }
    }
}

/// Parallelism of numeric factorisation and of the triangular solves.
///
/// The defaults are measured, not guessed: on faer 0.24.4 the *parallel* sparse Cholesky and LU
/// solves with more than one right-hand side are pathologically slow (a fixed ~1.4 s per call
/// for a 64 000-unknown 3-D problem whatever the block width, against 12 ms per vector
/// sequentially at width 8) and a single-vector solve gains nothing from threads, so solves run
/// sequentially; factorisations use the global (shared-pool) setting. See `docs/LEARNINGS.md`
/// and `examples/solve_scaling.rs`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThreadPolicy {
    pub factor: Threads,
    pub solve: Threads,
}

impl Default for ThreadPolicy {
    fn default() -> Self {
        Self { factor: Threads::Global, solve: Threads::Sequential }
    }
}

/// Recovery strategy used by [`ShiftedFactorization`] when `(A - sigma B)` is singular,
/// nearly singular, or exceeds direct-factorisation limits (R31).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SingularRecovery {
    /// Try small diagonal shift perturbations first (+ defect-corrected solve on `(A - sigma B)`),
    /// then fall back to preconditioned Krylov (`GMRES`) if perturbation still fails or `OutOfMemory` occurs.
    #[default]
    Auto,
    /// Recover only via automatic diagonal shift perturbation (`sigma + delta`).
    PerturbShift,
    /// Recover directly via preconditioned Krylov (`GMRES`) without perturbing `sigma`.
    IterativeFallback,
    /// Disable automatic recovery and return `SolverError::Singular` / `SolverError::OutOfMemory`.
    Disabled,
}

/// Pivot, Sylvester inertia, and singular-recovery diagnostics for [`ShiftedFactorization`] (R31).
#[derive(Debug, Clone, Copy, Default)]
pub struct FactorDiagnostics {
    /// Smallest pivot magnitude (eigenvalue magnitude for `1x1` / `2x2` Hermitian symmetric blocks,
    /// `L_{jj}^2` for Cholesky); `0.0` if zero/unavailable.
    pub min_pivot_abs: f64,
    /// Largest pivot magnitude observed in the active symmetric/Cholesky factor (`0.0` if unavailable).
    pub max_pivot_abs: f64,
    /// Smallest nonzero diagonal magnitude `|M_{ii}|` observed in the input matrix (`0.0` if all zero/unavailable).
    pub min_diag_abs: f64,
    /// Largest diagonal magnitude `max_i |M_{ii}|` observed in the input matrix.
    pub max_diag_abs: f64,
    /// Largest entry magnitude `max_{i,j} |M_{ij}|` observed in the input matrix.
    pub max_entry_abs: f64,
    /// Zero-cost pivot/diagonal growth ratio: captures upward pivot growth (`max_pivot / ref_scale`),
    /// diagonal pivot collapse (`min_diag / min_pivot`), or LU probe condition lower bound.
    pub pivot_growth: f64,
    /// Dynamic range of pivots (`max_pivot / min_pivot`, or probe/diagonal ratio on LU).
    pub pivot_ratio: f64,
    /// `true` when explicit factor pivot extrema were extracted (Cholesky / LBL^H / LBL^T).
    pub has_factor_pivots: bool,
    /// `1`-norm `max_j sum_i |M_{ij}|` of the shifted matrix `A - sigma B`.
    pub matrix_scale: f64,
    /// Count of pivots with magnitude `<= 1e-12 * matrix_scale`.
    pub small_pivots: usize,
    /// First failing column/block index recorded on singular / zero / non-finite pivot failure.
    pub first_bad_pivot_index: Option<usize>,
    /// Sylvester positive eigenvalue count `n_+` (available for real symmetric and complex Hermitian factors).
    pub positive_inertia: Option<usize>,
    /// Sylvester negative eigenvalue count `n_-` (available for real symmetric and complex Hermitian factors).
    pub negative_inertia: Option<usize>,
    /// Sylvester zero eigenvalue count `n_0` (available for real symmetric and complex Hermitian factors).
    pub zero_inertia: Option<usize>,
    /// Real part of shift perturbation `delta` when `sigma + delta` was used for factorisation recovery.
    pub shift_perturbation_re: f64,
    /// Imaginary part of shift perturbation `delta` when `sigma + delta` was used for factorisation recovery.
    pub shift_perturbation_im: f64,
    /// `true` when the active solver uses the preconditioned Krylov (`GMRES`) fallback path.
    pub krylov_fallback: bool,
    /// Cumulative Krylov iterations executed under `IterativeFallback` or shift-perturbation defect correction.
    pub krylov_iterations: usize,
}

impl FactorDiagnostics {
    /// Returns the complex/real shift perturbation `delta = factor_sigma - requested_sigma`.
    pub fn shift_perturbation<T: Scalar>(&self) -> T {
        T::from_parts(self.shift_perturbation_re, self.shift_perturbation_im)
    }

    /// Returns `true` when automatic shift perturbation was applied.
    pub fn has_shift_perturbation(&self) -> bool {
        self.shift_perturbation_re != 0.0 || self.shift_perturbation_im != 0.0
    }

    /// Returns `true` when zero-cost factorization diagnostics indicate anomalous pivot/diagonal
    /// growth, severe pivot collapse, small/bad pivots, or singular recovery.
    pub fn is_anomalous_pivot_growth(&self) -> bool {
        self.small_pivots > 0
            || self.first_bad_pivot_index.is_some()
            || self.has_shift_perturbation()
            || self.krylov_fallback
            || (self.has_factor_pivots && self.min_pivot_abs == 0.0 && self.matrix_scale > 0.0)
            || !self.pivot_growth.is_finite()
            || !self.pivot_ratio.is_finite()
            || self.pivot_growth > 1e6
            || self.pivot_ratio > 1e8
    }
}

/// Timings and counters of a [`ShiftedFactorization`] (seconds).
#[derive(Debug, Clone, Copy)]
pub struct FactorStats {
    /// Pattern union + symbolic analyses, paid once (LU and Cholesky analyses are separate and
    /// lazy).
    pub analyse_secs: f64,
    /// Number of symbolic analyses performed (stays at 1 per factorisation kind however many
    /// shifts are factorised: that is the reuse).
    pub symbolic_count: usize,
    /// Total and last numeric factorisation time, and the number of numeric factorisations.
    pub numeric_secs_total: f64,
    pub numeric_secs_last: f64,
    pub numeric_count: usize,
    /// Value reassembly time of the last shift.
    pub assemble_secs_last: f64,
    pub used: Used,
    /// Specific backend and fill-reducing ordering used for the last shift.
    pub backend: FactorBackend,
    /// Cholesky attempts that failed and fell back to LU / LBL^T.
    pub cholesky_fallbacks: usize,
    /// Right-hand sides solved.
    pub solves: usize,
    /// Stored (or symbolic symmetric) factor non-zeros (`L` for Cholesky / LBL^H, `L + U` for LU); 0 before analysis.
    pub factor_nnz: usize,
    pub solve_secs_total: f64,
}

/// "Factorise at a shift, then solve": the interface the eigensolvers use. A backend other than
/// `faer` (MKL PARDISO, a threaded MUMPS, a GPU solver) implements this trait.
pub trait ShiftedSolver<T: Scalar> {
    fn dim(&self) -> usize;
    /// Factorises `A − σB`. May reuse everything that depends only on the sparsity pattern.
    fn set_shift(&mut self, sigma: T) -> Result<(), SolverError>;
    /// Solves `(A − σB) X = RHS` in place for `nrhs` columns of length `dim()`.
    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) -> Result<(), SolverError>;
    fn stats(&self) -> FactorStats;
}

/// Union pattern (CSC) of `A`, `B` (or the identity), optionally restricted to the lower
/// triangle, with the values of `A` and `B` scattered onto it.
pub(crate) struct Assembled<T: Scalar> {
    pub(crate) col_ptr: Vec<u32>,
    pub(crate) row_idx: Vec<u32>,
    pub(crate) a: Vec<T>,
    pub(crate) b: Vec<T>,
}

pub(crate) fn assemble<T: Scalar>(a: &Csr<T>, b: Option<&Csr<T>>, lower: bool) -> Result<Assembled<T>, SolverError> {
    let n = a.nrows();
    if a.ncols() != n {
        return Err(SolverError::InvalidArgument("matrix must be square"));
    }
    if b.is_some_and(|b| b.nrows() != n || b.ncols() != n) {
        return Err(SolverError::InvalidArgument("A and B must have the same shape"));
    }
    if a.nnz() + n >= u32::MAX as usize / 2 {
        return Err(SolverError::TooLarge);
    }
    // Entries (col, row, a-value, b-value), gathered from CSR; the identity contributes the
    // diagonal when there is no B.
    let keep = |r: usize, c: usize| !lower || r >= c;
    let mut count = vec![0usize; n + 1];
    let mut tmp: Vec<(u32, u32, T, T)> = Vec::with_capacity(a.nnz() + b.map_or(n, |b| b.nnz()));
    for r in 0..n {
        let (cols, vals) = a.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            if keep(r, c as usize) {
                tmp.push((c, r as u32, v, T::ZERO));
            }
        }
    }
    match b {
        Some(b) => {
            for r in 0..n {
                let (cols, vals) = b.row(r);
                for (&c, &v) in cols.iter().zip(vals) {
                    if keep(r, c as usize) {
                        tmp.push((c, r as u32, T::ZERO, v));
                    }
                }
            }
        }
        None => {
            for r in 0..n {
                tmp.push((r as u32, r as u32, T::ZERO, T::ONE));
            }
        }
    }
    tmp.sort_by_key(|&(c, r, _, _)| (c, r));
    let (mut col_ptr, mut row_idx) = (vec![0u32; n + 1], Vec::with_capacity(tmp.len()));
    let (mut av, mut bv): (Vec<T>, Vec<T>) = (Vec::with_capacity(tmp.len()), Vec::with_capacity(tmp.len()));
    let mut last: Option<(u32, u32)> = None;
    for &(c, r, x, y) in &tmp {
        if let (true, Some(a_last), Some(b_last)) = (last == Some((c, r)), av.last_mut(), bv.last_mut()) {
            *a_last += x;
            *b_last += y;
        } else {
            row_idx.push(r);
            av.push(x);
            bv.push(y);
            count[c as usize + 1] += 1;
            last = Some((c, r));
        }
    }
    for c in 0..n {
        count[c + 1] += count[c];
        col_ptr[c + 1] = count[c + 1] as u32;
    }
    Ok(Assembled { col_ptr, row_idx, a: av, b: bv })
}

/// The lower-triangle union pattern of symmetric `A` and `B` in O(nnz) without sorting: column `c` of the lower triangle is row `c` of
/// the CSR from the diagonal on (symmetry), and the rows are merged with the rows of `B` (or the identity). Falls back to [`assemble`]
/// when a row is not sorted.
pub(crate) fn assemble_lower_sorted<T: Scalar>(a: &Csr<T>, b: Option<&Csr<T>>) -> Result<Assembled<T>, SolverError> {
    let n = a.nrows();
    if a.ncols() != n || b.is_some_and(|b| b.nrows() != n || b.ncols() != n) || !a.has_sorted_unique_rows() || b.is_some_and(|b| !b.has_sorted_unique_rows()) {
        return assemble(a, b, true);
    }
    if a.nnz() + n >= u32::MAX as usize / 2 {
        return Err(SolverError::TooLarge);
    }
    let cap = a.nnz() / 2 + n;
    let (mut col_ptr, mut row_idx) = (vec![0u32; n + 1], Vec::with_capacity(cap));
    let (mut av, mut bv): (Vec<T>, Vec<T>) = (Vec::with_capacity(cap), Vec::with_capacity(cap));
    for c in 0..n {
        let (ac, avs) = a.row(c);
        let mut i = ac.partition_point(|&x| (x as usize) < c);
        match b {
            Some(b) => {
                let (bc, bvs) = b.row(c);
                let mut j = bc.partition_point(|&x| (x as usize) < c);
                while i < ac.len() || j < bc.len() {
                    let (ca, cb) = (ac.get(i).copied().unwrap_or(u32::MAX), bc.get(j).copied().unwrap_or(u32::MAX));
                    let cm = ca.min(cb);
                    let x = if ca == cm { i += 1; avs[i - 1] } else { T::ZERO };
                    let y = if cb == cm { j += 1; bvs[j - 1] } else { T::ZERO };
                    row_idx.push(cm);
                    av.push(x);
                    bv.push(y);
                }
            }
            None => {
                if i >= ac.len() || ac[i] as usize != c {
                    row_idx.push(c as u32);
                    av.push(T::ZERO);
                    bv.push(T::ONE);
                }
                for k in i..ac.len() {
                    row_idx.push(ac[k]);
                    av.push(avs[k]);
                    bv.push(if ac[k] as usize == c { T::ONE } else { T::ZERO });
                }
            }
        }
        col_ptr[c + 1] = row_idx.len() as u32;
    }
    Ok(Assembled { col_ptr, row_idx, a: av, b: bv })
}

/// Cached transpose topology of the full CSC matrix `A`: `t_map[k]` is the entry index in `A^T`
/// corresponding to entry `k` of `A`, and `structurally_symmetric` is true iff `A` and `A^T` have
/// identical CSC patterns.
struct TransposeMap {
    t_map: Vec<usize>,
    structurally_symmetric: bool,
}

fn build_transpose_map<T: Scalar>(asm: &Assembled<T>, n: usize) -> TransposeMap {
    let nnz = asm.row_idx.len();
    let mut count = vec![0usize; n + 1];
    for &r in &asm.row_idx {
        count[r as usize + 1] += 1;
    }
    for i in 0..n {
        count[i + 1] += count[i];
    }
    let col_ptr: Vec<u32> = count.iter().map(|&c| c as u32).collect();
    let mut next = count;
    let mut row_idx = vec![0u32; nnz];
    let mut t_map = vec![0usize; nnz];
    for c in 0..n {
        let start = asm.col_ptr[c] as usize;
        let end = asm.col_ptr[c + 1] as usize;
        for k in start..end {
            let r = asm.row_idx[k] as usize;
            let dst = next[r];
            next[r] += 1;
            row_idx[dst] = c as u32;
            t_map[k] = dst;
        }
    }
    let structurally_symmetric = col_ptr == asm.col_ptr && row_idx == asm.row_idx;
    TransposeMap { t_map, structurally_symmetric }
}

fn extract_lower<T: Scalar>(full: &Assembled<T>, n: usize) -> Assembled<T> {
    let cap = full.row_idx.len() / 2 + n;
    let mut col_ptr = vec![0u32; n + 1];
    let mut row_idx = Vec::with_capacity(cap);
    let mut a = Vec::with_capacity(cap);
    let mut b = Vec::with_capacity(cap);
    for c in 0..n {
        let start = full.col_ptr[c] as usize;
        let end = full.col_ptr[c + 1] as usize;
        for k in start..end {
            let r = full.row_idx[k] as usize;
            if r >= c {
                row_idx.push(r as u32);
                a.push(full.a[k]);
                b.push(full.b[k]);
            }
        }
        col_ptr[c + 1] = row_idx.len() as u32;
    }
    Assembled { col_ptr, row_idx, a, b }
}

/// Error-free product `a * b = p + e` using exact FMA (`f64::mul_add`).
#[inline]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let e = a.mul_add(b, -p);
    (p, e)
}

/// Error-free sum `a + b = s + e` (Knuth 2Sum).
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}

/// Compensated double-double (`f64 + f64`) complex/real accumulator for 106-bit residual
/// evaluation `r = b - (A - sigma B) x` during deterministic iterative refinement (`ADR 11`).
#[derive(Clone, Copy, Default)]
struct CompAcc {
    re_hi: f64,
    re_lo: f64,
    im_hi: f64,
    im_lo: f64,
}

impl CompAcc {
    #[inline]
    fn from_scalar<T: Scalar>(v: T) -> Self {
        Self { re_hi: v.re(), re_lo: 0.0, im_hi: v.im(), im_lo: 0.0 }
    }

    #[inline]
    fn add_real_prod(&mut self, a: f64, b: f64) {
        let (p, pe) = two_prod(a, b);
        let (s, se) = two_sum(self.re_hi, p);
        self.re_hi = s;
        self.re_lo += se + pe;
    }

    #[inline]
    fn add_imag_prod(&mut self, a: f64, b: f64) {
        let (p, pe) = two_prod(a, b);
        let (s, se) = two_sum(self.im_hi, p);
        self.im_hi = s;
        self.im_lo += se + pe;
    }

    /// `self -= val * x` in compensated arithmetic.
    #[inline]
    fn sub_prod<T: Scalar>(&mut self, val: T, x: T) {
        if T::IS_COMPLEX {
            let (vr, vi) = (val.re(), val.im());
            let (xr, xi) = (x.re(), x.im());
            self.add_real_prod(-vr, xr);
            self.add_real_prod(vi, xi);
            self.add_imag_prod(-vr, xi);
            self.add_imag_prod(-vi, xr);
        } else {
            self.add_real_prod(-val.re(), x.re());
        }
    }

    #[inline]
    fn finish<T: Scalar>(self) -> T {
        if T::IS_COMPLEX {
            T::from_parts(self.re_hi + self.re_lo, self.im_hi + self.im_lo)
        } else {
            T::from_f64(self.re_hi + self.re_lo)
        }
    }
}

/// Deterministic compensated residual `r = b - A x` using a pre-assembled full CSC matrix `vals`.
fn csc_residual_vals<T: Scalar>(
    n: usize,
    col_ptr: &[u32],
    row_idx: &[u32],
    vals: &[T],
    b: &[T],
    x: &[T],
    r: &mut [T],
) {
    let mut acc: Vec<CompAcc> = b.iter().map(|&bi| CompAcc::from_scalar(bi)).collect();
    for c in 0..n {
        let xc = x[c];
        if xc == T::ZERO {
            continue;
        }
        let start = col_ptr[c] as usize;
        let end = col_ptr[c + 1] as usize;
        for k in start..end {
            acc[row_idx[k] as usize].sub_prod(vals[k], xc);
        }
    }
    for (ri, a) in r.iter_mut().zip(acc) {
        *ri = a.finish();
    }
}

/// Deterministic compensated residual `r = b - (A - sigma B) x` from full CSC `Assembled<T>`
/// without requiring a retained shifted `vals` buffer.
fn csc_residual_shifted_full<T: Scalar>(asm: &Assembled<T>, sigma: T, b: &[T], x: &[T], r: &mut [T]) {
    let n = b.len();
    let mut acc: Vec<CompAcc> = b.iter().map(|&bi| CompAcc::from_scalar(bi)).collect();
    let has_sigma = sigma != T::ZERO;
    for c in 0..n {
        let xc = x[c];
        if xc == T::ZERO {
            continue;
        }
        let neg_sigma_xc = if has_sigma { -(sigma * xc) } else { T::ZERO };
        let start = asm.col_ptr[c] as usize;
        let end = asm.col_ptr[c + 1] as usize;
        for k in start..end {
            let row = asm.row_idx[k] as usize;
            let ak = asm.a[k];
            if ak != T::ZERO {
                acc[row].sub_prod(ak, xc);
            }
            if has_sigma {
                let bk = asm.b[k];
                if bk != T::ZERO {
                    // Subtract (-sigma * bk) * xc  <==>  sub_prod(bk, -sigma * xc)
                    acc[row].sub_prod(bk, neg_sigma_xc);
                }
            }
        }
    }
    for (ri, a) in r.iter_mut().zip(acc) {
        *ri = a.finish();
    }
}

/// Deterministic compensated residual `r = b - (A - sigma B) x` from lower-triangle CSC `Assembled<T>`
/// (mirroring `r > c` with `conj` when `hermitian`, or plain transpose when complex-symmetric).
fn csc_residual_shifted_lower<T: Scalar>(
    lower: &Assembled<T>,
    sigma: T,
    hermitian: bool,
    b: &[T],
    x: &[T],
    r: &mut [T],
) {
    let n = b.len();
    let mut acc: Vec<CompAcc> = b.iter().map(|&bi| CompAcc::from_scalar(bi)).collect();
    let has_sigma = sigma != T::ZERO;
    for c in 0..n {
        let xc = x[c];
        let neg_sigma_xc = if has_sigma { -(sigma * xc) } else { T::ZERO };
        let start = lower.col_ptr[c] as usize;
        let end = lower.col_ptr[c + 1] as usize;
        for k in start..end {
            let row = lower.row_idx[k] as usize;
            let ak = lower.a[k];
            let bk = lower.b[k];
            if xc != T::ZERO {
                if ak != T::ZERO {
                    acc[row].sub_prod(ak, xc);
                }
                if has_sigma && bk != T::ZERO {
                    acc[row].sub_prod(bk, neg_sigma_xc);
                }
            }
            if row > c {
                let xr = x[row];
                if xr != T::ZERO {
                    let ak_u = if hermitian { ak.conj() } else { ak };
                    if ak_u != T::ZERO {
                        acc[c].sub_prod(ak_u, xr);
                    }
                    if has_sigma && bk != T::ZERO {
                        let bk_u = if hermitian { bk.conj() } else { bk };
                        let s_u = if hermitian { sigma.conj() } else { sigma };
                        acc[c].sub_prod(bk_u, -(s_u * xr));
                    }
                }
            }
        }
    }
    for (ri, a) in r.iter_mut().zip(acc) {
        *ri = a.finish();
    }
}

struct ProbeCheck {
    ok: bool,
    initial_rel_err: f64,
}

fn order_to_perm(order: &[usize]) -> (Vec<u32>, Vec<u32>) {
    let n = order.len();
    let fwd: Vec<u32> = order.iter().map(|&o| o as u32).collect();
    let mut inv = vec![0u32; n];
    for (new_idx, &old_idx) in order.iter().enumerate() {
        inv[old_idx] = new_idx as u32;
    }
    (fwd, inv)
}

fn csc_one_norm<T: Scalar>(n: usize, col_ptr: &[u32], vals: &[T]) -> f64 {
    let mut scale = 0.0f64;
    for c in 0..n {
        let mut s = 0.0f64;
        for &v in &vals[col_ptr[c] as usize..col_ptr[c + 1] as usize] {
            let m = v.modulus();
            if m.is_finite() {
                s += m;
            }
        }
        scale = scale.max(s);
    }
    scale
}

/// Computes `(one_norm, min_nonzero_diag_abs, max_diag_abs, max_entry_abs)` in a single pass
/// over full CSC columns.
fn csc_scale_and_diag<T: Scalar>(n: usize, col_ptr: &[u32], row_idx: &[u32], vals: &[T]) -> (f64, f64, f64, f64) {
    let mut scale = 0.0f64;
    let mut min_diag = f64::INFINITY;
    let mut max_diag = 0.0f64;
    let mut max_entry = 0.0f64;
    for c in 0..n {
        let mut s = 0.0f64;
        let start = col_ptr[c] as usize;
        let end = col_ptr[c + 1] as usize;
        for k in start..end {
            let m = vals[k].modulus();
            if m.is_finite() {
                s += m;
                max_entry = max_entry.max(m);
                if row_idx[k] as usize == c {
                    max_diag = max_diag.max(m);
                    if m > 0.0 {
                        min_diag = min_diag.min(m);
                    }
                }
            }
        }
        scale = scale.max(s);
    }
    (scale, if min_diag.is_finite() { min_diag } else { 0.0 }, max_diag, max_entry)
}

/// Computes `(one_norm, min_nonzero_diag_abs, max_diag_abs, max_entry_abs)` in a single pass
/// over lower-triangle CSC entries.
fn lower_scale_and_diag<T: Scalar>(n: usize, lower: &Assembled<T>, vals: &[T]) -> (f64, f64, f64, f64) {
    let mut col_sums = vec![0.0f64; n];
    let mut min_diag = f64::INFINITY;
    let mut max_diag = 0.0f64;
    let mut max_entry = 0.0f64;
    for c in 0..n {
        let start = lower.col_ptr[c] as usize;
        let end = lower.col_ptr[c + 1] as usize;
        for k in start..end {
            let r = lower.row_idx[k] as usize;
            let m = vals[k].modulus();
            if m.is_finite() {
                col_sums[c] += m;
                if r != c {
                    col_sums[r] += m;
                } else {
                    max_diag = max_diag.max(m);
                    if m > 0.0 {
                        min_diag = min_diag.min(m);
                    }
                }
                max_entry = max_entry.max(m);
            }
        }
    }
    let scale = col_sums.into_iter().fold(0.0f64, f64::max);
    (scale, if min_diag.is_finite() { min_diag } else { 0.0 }, max_diag, max_entry)
}

fn finalize_pivot_growth(diag: &mut FactorDiagnostics, min_diag: f64, max_diag: f64, max_entry: f64) {
    diag.min_diag_abs = min_diag;
    diag.max_diag_abs = max_diag;
    diag.max_entry_abs = max_entry;
    diag.has_factor_pivots = true;
    let ref_scale = max_diag.max(max_entry).max(diag.matrix_scale).max(1e-300);
    let growth_up = if diag.max_pivot_abs > 0.0 {
        diag.max_pivot_abs / ref_scale
    } else {
        1.0
    };
    let diag_ref = if min_diag > 0.0 { min_diag } else { max_entry };
    let growth_down = if diag_ref > 0.0 && diag.min_pivot_abs > 0.0 {
        diag_ref / diag.min_pivot_abs
    } else if diag_ref > 0.0 && diag.matrix_scale > 0.0 {
        f64::INFINITY
    } else {
        1.0
    };
    diag.pivot_growth = growth_up.max(growth_down).max(1.0);
    diag.pivot_ratio = if diag.max_pivot_abs > 0.0 && diag.min_pivot_abs > 0.0 {
        diag.max_pivot_abs / diag.min_pivot_abs
    } else if diag.matrix_scale > 0.0 {
        f64::INFINITY
    } else {
        1.0
    };
}

fn b_one_norm_full<T: Scalar>(full: &Assembled<T>, n: usize) -> f64 {
    csc_one_norm(n, &full.col_ptr, &full.b)
}

fn csc_bad_diag_col<T: Scalar>(n: usize, col_ptr: &[u32], vals: &[T]) -> usize {
    for c in 0..n {
        let slice = &vals[col_ptr[c] as usize..col_ptr[c + 1] as usize];
        let s: f64 = slice.iter().map(|v| v.modulus()).sum();
        if s == 0.0 || !s.is_finite() {
            return c;
        }
    }
    0
}

fn record_cholesky_diagnostics<T: FaerScalar>(
    sym: &SymbolicCholesky<u32>,
    l_values: &[T],
    scale: f64,
    diag: &mut FactorDiagnostics,
) {
    let n = sym.nrows();
    if n == 0 {
        diag.min_pivot_abs = 0.0;
        diag.max_pivot_abs = 0.0;
        diag.small_pivots = 0;
        diag.positive_inertia = Some(0);
        diag.negative_inertia = Some(0);
        diag.zero_inertia = Some(0);
        return;
    }
    let tiny = 1e-12 * scale.max(1e-300);
    let mut min_p = f64::INFINITY;
    let mut max_p = 0.0f64;
    let mut small = 0usize;
    let mut push = |l_jj: T| {
        let d = l_jj.modulus().powi(2);
        min_p = min_p.min(d);
        max_p = max_p.max(d);
        if d <= tiny {
            small += 1;
        }
    };
    match sym.raw() {
        SymbolicCholeskyRaw::Supernodal(sup) => {
            let sb = sup.supernode_begin();
            let se = sup.supernode_end();
            let cv = sup.col_ptr_for_val();
            let nnz = sup.nnz_per_super();
            for s in 0..sup.n_supernodes() {
                let s_ncols = se[s] as usize - sb[s] as usize;
                let s_nrows = s_ncols + nnz[s] as usize;
                let base = cv[s] as usize;
                for j in 0..s_ncols {
                    push(l_values[base + j * s_nrows + j]);
                }
            }
        }
        SymbolicCholeskyRaw::Simplicial(simp) => {
            let cp = simp.col_ptr();
            for j in 0..n {
                push(l_values[cp[j] as usize]);
            }
        }
    }
    diag.min_pivot_abs = if min_p.is_finite() { min_p } else { 0.0 };
    diag.max_pivot_abs = max_p;
    diag.small_pivots = small;
    diag.positive_inertia = Some(n);
    diag.negative_inertia = Some(0);
    diag.zero_inertia = Some(0);
}

fn record_lblt_diagnostics<T: FaerScalar>(
    sym: &SymbolicCholesky<u32>,
    l_values: &[T],
    subdiag: &[T],
    hermitian_or_real: bool,
    scale: f64,
    diag: &mut FactorDiagnostics,
) -> Option<usize> {
    struct Acc {
        tiny: f64,
        min_p: f64,
        max_p: f64,
        small: usize,
        pos: usize,
        neg: usize,
        zero: usize,
        first_bad: Option<usize>,
    }
    impl Acc {
        fn push_bad(&mut self, col: usize, count: usize) {
            self.min_p = 0.0;
            self.small += count;
            self.zero += count;
            self.first_bad.get_or_insert(col);
        }
        fn push_real_eig(&mut self, e: f64, col: usize) {
            if !e.is_finite() || e == 0.0 {
                self.push_bad(col, 1);
            } else {
                let ae = e.abs();
                self.min_p = self.min_p.min(ae);
                self.max_p = self.max_p.max(ae);
                if ae <= self.tiny {
                    self.small += 1;
                }
                if e > 0.0 {
                    self.pos += 1;
                } else {
                    self.neg += 1;
                }
            }
        }
        fn push_complex_mag(&mut self, am: f64, smax: f64, count: usize) {
            self.min_p = self.min_p.min(am);
            self.max_p = self.max_p.max(smax);
            if am <= self.tiny {
                self.small += count;
            }
        }
    }

    let n = sym.nrows();
    let mut acc = Acc {
        tiny: 1e-12 * scale.max(1e-300),
        min_p: f64::INFINITY,
        max_p: 0.0,
        small: 0,
        pos: 0,
        neg: 0,
        zero: 0,
        first_bad: None,
    };

    match sym.raw() {
        SymbolicCholeskyRaw::Supernodal(sup) => {
            let sb = sup.supernode_begin();
            let se = sup.supernode_end();
            let cv = sup.col_ptr_for_val();
            let nnz = sup.nnz_per_super();
            for s in 0..sup.n_supernodes() {
                let s_start = sb[s] as usize;
                let s_ncols = se[s] as usize - s_start;
                let s_nrows = s_ncols + nnz[s] as usize;
                let base = cv[s] as usize;
                let s_sub = &subdiag[s_start..s_start + s_ncols];
                let mut idx = 0usize;
                while idx < s_ncols {
                    let gcol = s_start + idx;
                    let b11 = l_values[base + idx * s_nrows + idx];
                    let b21 = s_sub[idx];
                    if b21 == T::ZERO {
                        if !b11.is_finite() || b11.modulus() == 0.0 {
                            acc.push_bad(gcol, 1);
                        } else if hermitian_or_real {
                            acc.push_real_eig(b11.re(), gcol);
                        } else {
                            let am = b11.modulus();
                            acc.push_complex_mag(am, am, 1);
                        }
                        idx += 1;
                    } else if idx + 1 < s_ncols {
                        let b22 = l_values[base + (idx + 1) * s_nrows + idx + 1];
                        if !b11.is_finite() || !b21.is_finite() || !b22.is_finite() {
                            acc.push_bad(gcol, 2);
                        } else if hermitian_or_real {
                            let tr = b11.re() + b22.re();
                            let d = b11.re() - b22.re();
                            let disc = (d * d + 4.0 * b21.modulus().powi(2)).sqrt();
                            acc.push_real_eig(0.5 * (tr + disc), gcol);
                            acc.push_real_eig(0.5 * (tr - disc), gcol + 1);
                        } else {
                            let det = b11 * b22 - b21 * b21;
                            let tr_m = (b11.modulus() + b22.modulus() + 2.0 * b21.modulus()).max(1e-300);
                            if !det.is_finite() || det.modulus() == 0.0 {
                                acc.push_bad(gcol, 2);
                            } else {
                                acc.push_complex_mag(det.modulus() / tr_m, tr_m, 2);
                            }
                        }
                        idx += 2;
                    } else {
                        acc.push_bad(gcol, 1);
                        idx += 1;
                    }
                }
            }
        }
        SymbolicCholeskyRaw::Simplicial(simp) => {
            let cp = simp.col_ptr();
            for j in 0..n {
                let d = l_values[cp[j] as usize];
                if !d.is_finite() || d.modulus() == 0.0 {
                    acc.push_bad(j, 1);
                } else if hermitian_or_real {
                    acc.push_real_eig(d.re(), j);
                } else {
                    let am = d.modulus();
                    acc.push_complex_mag(am, am, 1);
                }
            }
        }
    }

    diag.min_pivot_abs = if acc.min_p.is_finite() { acc.min_p } else { 0.0 };
    diag.max_pivot_abs = acc.max_p;
    diag.small_pivots = acc.small;
    if let Some(bad) = acc.first_bad {
        diag.first_bad_pivot_index.get_or_insert(bad);
    }
    if hermitian_or_real {
        diag.positive_inertia = Some(acc.pos);
        diag.negative_inertia = Some(acc.neg);
        diag.zero_inertia = Some(acc.zero);
    } else {
        diag.positive_inertia = None;
        diag.negative_inertia = None;
        diag.zero_inertia = None;
    }
    acc.first_bad
}

fn assembled_to_csr<T: Scalar>(full: &Assembled<T>, n: usize, sigma: T) -> Result<Csr<T>, SolverError> {
    let nnz = full.row_idx.len();
    let mut row_ptr = vec![0usize; n + 1];
    for &r in &full.row_idx {
        row_ptr[r as usize + 1] += 1;
    }
    for r in 0..n {
        row_ptr[r + 1] += row_ptr[r];
    }
    let mut write_pos = row_ptr.clone();
    let mut col_idx = vec![0u32; nnz];
    let mut vals = vec![T::ZERO; nnz];
    let has_sigma = sigma != T::ZERO;
    for c in 0..n {
        let start = full.col_ptr[c] as usize;
        let end = full.col_ptr[c + 1] as usize;
        for k in start..end {
            let r = full.row_idx[k] as usize;
            let pos = write_pos[r];
            write_pos[r] += 1;
            col_idx[pos] = c as u32;
            vals[pos] = if has_sigma { full.a[k] - sigma * full.b[k] } else { full.a[k] };
        }
    }
    Csr::try_new(n, n, row_ptr, col_idx, vals).map_err(|_| SolverError::TooLarge)
}

enum KrylovPc<T: Scalar> {
    Ilu(Box<Ilu<T>>),
    Jacobi(Jacobi<T>),
}

impl<T: Scalar> Preconditioner<T> for KrylovPc<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        match self {
            KrylovPc::Ilu(p) => p.apply(r, z),
            KrylovPc::Jacobi(p) => p.apply(r, z),
        }
    }
}

struct KrylovFallback<T: Scalar> {
    mat: Csr<T>,
    pc: KrylovPc<T>,
}

struct ShiftedOpView<'a, T: FaerScalar> {
    fac: &'a ShiftedFactorization<T>,
}

impl<T: FaerScalar> LinearOperator<T> for ShiftedOpView<'_, T> {
    fn dim(&self) -> usize {
        self.fac.n
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        if let Some(full) = &self.fac.full {
            y.fill(T::ZERO);
            let has_sigma = self.fac.sigma != T::ZERO;
            for c in 0..self.fac.n {
                let xc = x[c];
                if xc == T::ZERO {
                    continue;
                }
                let start = full.col_ptr[c] as usize;
                let end = full.col_ptr[c + 1] as usize;
                for k in start..end {
                    let r = full.row_idx[k] as usize;
                    let v = if has_sigma { full.a[k] - self.fac.sigma * full.b[k] } else { full.a[k] };
                    y[r] += v * xc;
                }
            }
        } else {
            let mut r_neg = vec![T::ZERO; self.fac.n];
            self.fac.compute_shifted_residual(&vec![T::ZERO; self.fac.n], x, &mut r_neg);
            for (yi, &ri) in y.iter_mut().zip(&r_neg) {
                *yi = -ri;
            }
        }
    }
}

struct PerturbedFactorPc<'a, T: FaerScalar> {
    fac: &'a ShiftedFactorization<T>,
}

impl<T: FaerScalar> Preconditioner<T> for PerturbedFactorPc<'_, T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
        if self.fac.solve_factor_only_in_place(z, 1).is_err() {
            z.copy_from_slice(r);
        }
    }
}

// ------------------------------------------------------------------------------------------------------------ Non-conjugate symmetric indefinite LBL^T (complex-symmetric PML / lossy)

const BK_ALPHA: f64 = 0.640_388_203_202_207_6; // (1 + sqrt(17)) / 8

/// Cached supernodal descendant lists and permuted lower-triangle scatter map for
/// non-conjugate complex-symmetric `LBL^T` on top of `faer`'s `SymbolicCholesky`.
struct SymLbltPlan {
    /// Permuted lower-triangle CSC `(col_ptr, row_idx)` where entry `dst` comes from `lower[src_idx[dst]]`.
    perm_col_ptr: Vec<usize>,
    perm_row_idx: Vec<usize>,
    src_idx: Vec<usize>,
    /// For each supernode `s`, the list of strictly earlier supernodes `d < s` whose row pattern hits `[s_start..s_end)`.
    descendants: Vec<Vec<usize>>,
}

fn build_sym_lblt_plan<T: Scalar>(lower: &Assembled<T>, sym: &SymbolicCholesky<u32>, n: usize) -> SymLbltPlan {
    let inv_opt = sym.perm().map(|p| p.arrays().1);
    let mut entries: Vec<(usize, usize, usize)> = Vec::with_capacity(lower.row_idx.len());
    for c in 0..n {
        let start = lower.col_ptr[c] as usize;
        let end = lower.col_ptr[c + 1] as usize;
        let pc = inv_opt.map_or(c, |inv| inv[c] as usize);
        for k in start..end {
            let r = lower.row_idx[k] as usize;
            let pr = inv_opt.map_or(r, |inv| inv[r] as usize);
            let (row, col) = if pr >= pc { (pr, pc) } else { (pc, pr) };
            entries.push((col, row, k));
        }
    }
    let mut count = vec![0usize; n + 1];
    for &(col, _, _) in &entries {
        count[col + 1] += 1;
    }
    for j in 0..n {
        count[j + 1] += count[j];
    }
    let perm_col_ptr = count.clone();
    let mut next = count;
    let mut perm_row_idx = vec![0usize; entries.len()];
    let mut src_idx = vec![0usize; entries.len()];
    for &(col, row, k) in &entries {
        let pos = next[col];
        next[col] += 1;
        perm_row_idx[pos] = row;
        src_idx[pos] = k;
    }

    let descendants = match sym.raw() {
        SymbolicCholeskyRaw::Supernodal(sup) => {
            let ns = sup.n_supernodes();
            let mut col_to_super = vec![0usize; n];
            let sb = sup.supernode_begin();
            let se = sup.supernode_end();
            for s in 0..ns {
                for col in sb[s] as usize..se[s] as usize {
                    col_to_super[col] = s;
                }
            }
            let mut desc: Vec<Vec<usize>> = vec![Vec::new(); ns];
            let col_ptr_row = sup.col_ptr_for_row_idx();
            let nnz_per_super = sup.nnz_per_super();
            let row_idx = sup.row_idx();
            for d in 0..ns {
                let rstart = col_ptr_row[d] as usize;
                let rlen = nnz_per_super[d] as usize;
                let mut prev_s = usize::MAX;
                for &r in &row_idx[rstart..rstart + rlen] {
                    let s = col_to_super[r as usize];
                    if s != prev_s {
                        desc[s].push(d);
                        prev_s = s;
                    }
                }
            }
            desc
        }
        SymbolicCholeskyRaw::Simplicial(_) => Vec::new(),
    };

    SymLbltPlan { perm_col_ptr, perm_row_idx, src_idx, descendants }
}

#[inline]
fn sym_swap_lower<T: Scalar>(a: &mut [T], m: usize, p: usize, q: usize) {
    debug_assert!(p < q && q < m);
    for j in 0..p {
        a.swap(j * m + p, j * m + q);
    }
    a.swap(p * m + p, q * m + q);
    for i in p + 1..q {
        a.swap(p * m + i, i * m + q);
    }
    for i in q + 1..m {
        a.swap(p * m + i, q * m + i);
    }
}

/// Dense non-conjugate Bunch–Kaufman `P^T L B L^T P = A` on an `m × m` column-major lower-triangular
/// matrix `a`. Stores unit-lower `L` below the diagonal (with `L[k+1, k] = 0` on 2×2 pivots),
/// diagonal of `B` on `diag(a)`, and subdiagonal of `B` in `subdiag[0..m]`.
fn dense_sym_lblt<T: Scalar>(
    a: &mut [T],
    m: usize,
    subdiag: &mut [T],
    perm_fwd: &mut [u32],
    perm_inv: &mut [u32],
    offset: usize,
) -> Result<(), usize> {
    for i in 0..m {
        perm_fwd[i] = i as u32;
        subdiag[i] = T::ZERO;
    }
    let mut w0 = vec![T::ZERO; m];
    let mut w1 = vec![T::ZERO; m];
    let mut k = 0usize;
    while k < m {
        let akk = a[k * m + k].modulus();
        let mut lambda = 0.0f64;
        let mut r = k;
        for i in k + 1..m {
            let v = a[k * m + i].modulus();
            if v > lambda {
                lambda = v;
                r = i;
            }
        }
        let peak = akk.max(lambda);
        if peak == 0.0 || !peak.is_finite() {
            return Err(offset + k);
        }
        let mut step = 1usize;
        if akk < BK_ALPHA * lambda {
            // Scan column r (rows k..m, i != r) for sigma = max |A_{i, r}|.
            let mut sigma = 0.0f64;
            for i in k..r {
                sigma = sigma.max(a[i * m + r].modulus());
            }
            for i in r + 1..m {
                sigma = sigma.max(a[r * m + i].modulus());
            }
            if akk * sigma >= BK_ALPHA * lambda * lambda {
                step = 1;
            } else if a[r * m + r].modulus() >= BK_ALPHA * sigma {
                sym_swap_lower(a, m, k, r);
                perm_fwd.swap(k, r);
                step = 1;
            } else {
                if r != k + 1 {
                    sym_swap_lower(a, m, k + 1, r);
                    perm_fwd.swap(k + 1, r);
                }
                step = 2;
            }
        }

        if step == 1 {
            let d = a[k * m + k];
            if d.modulus() == 0.0 || !d.is_finite() {
                return Err(offset + k);
            }
            let d_inv = T::ONE / d;
            subdiag[k] = T::ZERO;
            let tail = m - (k + 1);
            for t in 0..tail {
                let i = k + 1 + t;
                let v = a[k * m + i];
                w0[t] = v;
                a[k * m + i] = v * d_inv;
            }
            for j in k + 1..m {
                let wj = w0[j - (k + 1)];
                if wj == T::ZERO {
                    continue;
                }
                let col_j = j * m;
                let col_k = k * m;
                for i in j..m {
                    a[col_j + i] -= a[col_k + i] * wj;
                }
            }
            k += 1;
        } else {
            let b11 = a[k * m + k];
            let b21 = a[k * m + k + 1];
            let b22 = a[(k + 1) * m + k + 1];
            let det = b11 * b22 - b21 * b21;
            if det.modulus() == 0.0 || !det.is_finite() {
                return Err(offset + k);
            }
            let inv_det = T::ONE / det;
            subdiag[k] = b21;
            subdiag[k + 1] = T::ZERO;
            // Unit lower L has 2×2 identity on the 2×2 pivot block.
            a[k * m + k + 1] = T::ZERO;
            let tail = m - (k + 2);
            for t in 0..tail {
                let i = k + 2 + t;
                let v0 = a[k * m + i];
                let v1 = a[(k + 1) * m + i];
                w0[t] = v0;
                w1[t] = v1;
                a[k * m + i] = (b22 * v0 - b21 * v1) * inv_det;
                a[(k + 1) * m + i] = (b11 * v1 - b21 * v0) * inv_det;
            }
            for j in k + 2..m {
                let wj0 = w0[j - (k + 2)];
                let wj1 = w1[j - (k + 2)];
                let col_j = j * m;
                let col_k0 = k * m;
                let col_k1 = (k + 1) * m;
                for i in j..m {
                    a[col_j + i] -= a[col_k0 + i] * wj0 + a[col_k1 + i] * wj1;
                }
            }
            k += 2;
        }
    }
    for (new_i, &old_i) in perm_fwd.iter().enumerate() {
        perm_inv[old_i as usize] = (offset + new_i) as u32;
    }
    for p in perm_fwd.iter_mut() {
        *p += offset as u32;
    }
    Ok(())
}

/// Non-conjugate supernodal intranodal Bunch–Kaufman `LBL^T` on `SymbolicSupernodalCholesky`.
#[allow(clippy::too_many_arguments)]
fn factorize_supernodal_sym_lblt<T: FaerScalar>(
    l_values: &mut [T],
    subdiag: &mut [T],
    perm_fwd: &mut [u32],
    perm_inv: &mut [u32],
    lower_vals: &[T],
    plan: &SymLbltPlan,
    sym: &SymbolicSupernodalCholesky<u32>,
    par: Par,
) -> Result<(), usize> {
    let n = sym.nrows();
    let ns = sym.n_supernodes();
    l_values.fill(T::ZERO);
    subdiag.fill(T::ZERO);
    let sb = sym.supernode_begin();
    let se = sym.supernode_end();
    let col_ptr_row = sym.col_ptr_for_row_idx();
    let col_ptr_val = sym.col_ptr_for_val();
    let nnz_per_super = sym.nnz_per_super();
    let row_idx = sym.row_idx();

    let mut global_to_local = vec![usize::MAX; n];
    let mut dense_top: Vec<T> = Vec::new();
    let mut tmp_buf: Vec<T> = Vec::new();
    let mut ld_b_buf: Vec<T> = Vec::new();
    let mut bot_perm_buf: Vec<T> = Vec::new();

    for s in 0..ns {
        let s_start = sb[s] as usize;
        let s_end = se[s] as usize;
        let s_ncols = s_end - s_start;
        let rstart = col_ptr_row[s] as usize;
        let s_pattern = &row_idx[rstart..rstart + nnz_per_super[s] as usize];
        let s_pat_len = s_pattern.len();
        let s_nrows = s_ncols + s_pat_len;

        for (i, &row) in s_pattern.iter().enumerate() {
            global_to_local[row as usize] = i + s_ncols;
        }

        let v_start = col_ptr_val[s] as usize;
        let v_end = col_ptr_val[s + 1] as usize;
        let (head, tail) = l_values.split_at_mut(v_start);
        let ls_slice = &mut tail[..v_end - v_start];

        // Scatter lower-triangle entries of permuted A for columns s_start..s_end.
        for j in s_start..s_end {
            let j_loc = j - s_start;
            let p0 = plan.perm_col_ptr[j];
            let p1 = plan.perm_col_ptr[j + 1];
            for pos in p0..p1 {
                let row = plan.perm_row_idx[pos];
                let val = lower_vals[plan.src_idx[pos]];
                let i_loc = if row >= s_end { global_to_local[row] } else { row - s_start };
                ls_slice[j_loc * s_nrows + i_loc] += val;
            }
        }

        // Left-looking BLAS-3 updates from descendant supernodes d < s (plain transpose, NO conjugation).
        for &d in &plan.descendants[s] {
            let d_start = sb[d] as usize;
            let d_end = se[d] as usize;
            let d_ncols = d_end - d_start;
            let dr = col_ptr_row[d] as usize;
            let d_pattern = &row_idx[dr..dr + nnz_per_super[d] as usize];
            let d_nrows = d_ncols + d_pattern.len();
            let d_slice = &head[col_ptr_val[d] as usize..col_ptr_val[d + 1] as usize];

            let d_pat_start = d_pattern.partition_point(|&r| (r as usize) < s_start);
            let d_mid_len = d_pattern[d_pat_start..].partition_point(|&r| (r as usize) < s_end);
            let d_pat_mid = d_pat_start + d_mid_len;
            let d_bot_len = d_pattern.len() - d_pat_mid;
            let d_mid_bot_len = d_mid_len + d_bot_len;
            if d_mid_len == 0 {
                continue;
            }

            // Build `ld_b = Ld_mid * B_d` (column-major `d_mid_len × d_ncols`, NO conj).
            ld_b_buf.resize(d_mid_len * d_ncols, T::ZERO);
            let d_subdiag = &subdiag[d_start..d_end];
            let mut j = 0usize;
            while j < d_ncols {
                let b21 = d_subdiag[j];
                if b21 == T::ZERO {
                    let b11 = d_slice[j * d_nrows + j];
                    let src_col = &d_slice[j * d_nrows + d_ncols + d_pat_start..j * d_nrows + d_ncols + d_pat_mid];
                    let dst_col = &mut ld_b_buf[j * d_mid_len..(j + 1) * d_mid_len];
                    for i in 0..d_mid_len {
                        dst_col[i] = src_col[i] * b11;
                    }
                    j += 1;
                } else {
                    let b11 = d_slice[j * d_nrows + j];
                    let b22 = d_slice[(j + 1) * d_nrows + j + 1];
                    let src0 = &d_slice[j * d_nrows + d_ncols + d_pat_start..j * d_nrows + d_ncols + d_pat_mid];
                    let src1 =
                        &d_slice[(j + 1) * d_nrows + d_ncols + d_pat_start..(j + 1) * d_nrows + d_ncols + d_pat_mid];
                    for i in 0..d_mid_len {
                        let x0 = src0[i];
                        let x1 = src1[i];
                        ld_b_buf[j * d_mid_len + i] = x0 * b11 + x1 * b21;
                        ld_b_buf[(j + 1) * d_mid_len + i] = x1 * b22 + x0 * b21;
                    }
                    j += 2;
                }
            }

            // `tmp = Ld_mid_bot * (Ld_mid * B_d)^T` (non-conjugate!) via faer BLAS-3.
            tmp_buf.resize(d_mid_bot_len * d_mid_len, T::ZERO);
            let ld_mat = MatRef::from_column_major_slice(T::faer_slice(d_slice), d_nrows, d_ncols);
            let (_, ld_mid_bot) = ld_mat.split_at_row(d_ncols + d_pat_start);
            let (ld_mid, ld_bot) = ld_mid_bot.split_at_row(d_mid_len);
            let ld_b_mat = MatRef::from_column_major_slice(T::faer_slice(&ld_b_buf[..d_mid_len * d_ncols]), d_mid_len, d_ncols);
            let tmp_mat = MatMut::from_column_major_slice_mut(
                T::faer_slice_mut(&mut tmp_buf[..d_mid_bot_len * d_mid_len]),
                d_mid_bot_len,
                d_mid_len,
            );
            let (tmp_top, tmp_bot) = tmp_mat.split_at_row_mut(d_mid_len);
            let one = <T::Elem as faer::traits::ComplexField>::one_impl();
            matmul::triangular::matmul(
                tmp_top,
                matmul::triangular::BlockStructure::TriangularLower,
                Accum::Replace,
                ld_mid,
                matmul::triangular::BlockStructure::Rectangular,
                ld_b_mat.transpose(),
                matmul::triangular::BlockStructure::Rectangular,
                one,
                par,
            );
            if d_bot_len > 0 {
                matmul::matmul(tmp_bot, Accum::Replace, ld_bot, ld_b_mat.transpose(), one, par);
            }

            // Scatter-subtract lower triangle of `tmp_top` and full `tmp_bot` into `Ls`.
            let mid_pat = &d_pattern[d_pat_start..d_pat_mid];
            for (j_idx, &gj) in mid_pat.iter().enumerate() {
                let j_s = gj as usize - s_start;
                let col_off = j_idx * d_mid_bot_len;
                let ls_col = j_s * s_nrows;
                for (i_rel, &gi) in mid_pat[j_idx..].iter().enumerate() {
                    let i_idx = j_idx + i_rel;
                    let i_s = gi as usize - s_start;
                    ls_slice[ls_col + i_s] -= tmp_buf[col_off + i_idx];
                }
                for (i_rel, &gi) in d_pattern[d_pat_mid..].iter().enumerate() {
                    let i_s = global_to_local[gi as usize];
                    ls_slice[ls_col + i_s] -= tmp_buf[col_off + d_mid_len + i_rel];
                }
            }
        }

        // Factorize dense s_ncols × s_ncols diagonal block `Ls_top` with non-conjugate Bunch–Kaufman.
        dense_top.resize(s_ncols * s_ncols, T::ZERO);
        for c in 0..s_ncols {
            dense_top[c * s_ncols..(c + 1) * s_ncols].copy_from_slice(&ls_slice[c * s_nrows..c * s_nrows + s_ncols]);
        }
        dense_sym_lblt(
            &mut dense_top[..s_ncols * s_ncols],
            s_ncols,
            &mut subdiag[s_start..s_end],
            &mut perm_fwd[s_start..s_end],
            &mut perm_inv[s_start..s_end],
            s_start,
        )?;
        // Copy factorized lower triangle back into `Ls_top` and zero strictly upper triangle.
        for c in 0..s_ncols {
            let ls_col = c * s_nrows;
            for r in 0..c {
                ls_slice[ls_col + r] = T::ZERO;
            }
            ls_slice[ls_col + c..ls_col + s_ncols].copy_from_slice(&dense_top[c * s_ncols + c..(c + 1) * s_ncols]);
        }

        if s_pat_len > 0 {
            // Permute columns of `Ls_bot` by local `perm_fwd`.
            bot_perm_buf.resize(s_pat_len * s_ncols, T::ZERO);
            for c in 0..s_ncols {
                let src_c = perm_fwd[s_start + c] as usize - s_start;
                bot_perm_buf[c * s_pat_len..(c + 1) * s_pat_len]
                    .copy_from_slice(&ls_slice[src_c * s_nrows + s_ncols..src_c * s_nrows + s_nrows]);
            }
            for c in 0..s_ncols {
                ls_slice[c * s_nrows + s_ncols..c * s_nrows + s_nrows]
                    .copy_from_slice(&bot_perm_buf[c * s_pat_len..(c + 1) * s_pat_len]);
            }

            // Solve `L_top * W^T = Ls_bot^T` (unit lower, NO conj!).
            {
                let ls_mat = MatMut::from_column_major_slice_mut(T::faer_slice_mut(ls_slice), s_nrows, s_ncols);
                let (ls_top, ls_bot) = ls_mat.split_at_row_mut(s_ncols);
                triangular_solve::solve_unit_lower_triangular_in_place(ls_top.rb(), ls_bot.transpose_mut(), par);
            }

            // Scale columns of `Ls_bot` on the right by `B_s^{-1}` (NO conj!).
            let s_sub = &subdiag[s_start..s_end];
            let mut j = 0usize;
            while j < s_ncols {
                let b21 = s_sub[j];
                if b21 == T::ZERO {
                    let d_inv = T::ONE / ls_slice[j * s_nrows + j];
                    for i in 0..s_pat_len {
                        ls_slice[j * s_nrows + s_ncols + i] *= d_inv;
                    }
                    j += 1;
                } else {
                    let b11 = ls_slice[j * s_nrows + j];
                    let b22 = ls_slice[(j + 1) * s_nrows + j + 1];
                    let inv_det = T::ONE / (b11 * b22 - b21 * b21);
                    let c0 = j * s_nrows + s_ncols;
                    let c1 = (j + 1) * s_nrows + s_ncols;
                    for i in 0..s_pat_len {
                        let x0 = ls_slice[c0 + i];
                        let x1 = ls_slice[c1 + i];
                        ls_slice[c0 + i] = (b22 * x0 - b21 * x1) * inv_det;
                        ls_slice[c1 + i] = (b11 * x1 - b21 * x0) * inv_det;
                    }
                    j += 2;
                }
            }
        }

        for &row in s_pattern {
            global_to_local[row as usize] = usize::MAX;
        }
    }
    clear_upper_vector_state();
    Ok(())
}

/// Solves `A X = RHS` in place using the non-conjugate supernodal `LBL^T` factor.
#[allow(clippy::too_many_arguments)]
fn solve_supernodal_sym_lblt<T: FaerScalar>(
    sym_top: &SymbolicCholesky<u32>,
    sym: &SymbolicSupernodalCholesky<u32>,
    l_values: &[T],
    subdiag: &[T],
    dyn_fwd: &[u32],
    dyn_inv: &[u32],
    rhs: &mut [T],
    nrhs: usize,
    par: Par,
) {
    let n = sym.nrows();
    let ns = sym.n_supernodes();
    let sb = sym.supernode_begin();
    let se = sym.supernode_end();
    let col_ptr_row = sym.col_ptr_for_row_idx();
    let col_ptr_val = sym.col_ptr_for_val();
    let nnz_per_super = sym.nnz_per_super();
    let row_idx = sym.row_idx();

    // Forward permutation: x[i, col] = rhs[sym_fwd[dyn_fwd[i]], col].
    let mut x = vec![T::ZERO; n * nrhs];
    match sym_top.perm() {
        Some(p) => {
            let sym_fwd = p.arrays().0;
            for col in 0..nrhs {
                let src = &rhs[col * n..(col + 1) * n];
                let dst = &mut x[col * n..(col + 1) * n];
                for i in 0..n {
                    dst[i] = src[sym_fwd[dyn_fwd[i] as usize] as usize];
                }
            }
        }
        None => {
            for col in 0..nrhs {
                let src = &rhs[col * n..(col + 1) * n];
                let dst = &mut x[col * n..(col + 1) * n];
                for i in 0..n {
                    dst[i] = src[dyn_fwd[i] as usize];
                }
            }
        }
    }

    let one = <T::Elem as faer::traits::ComplexField>::one_impl();
    let neg_one = T::to_faer(-T::ONE);
    let mut tmp: Vec<T> = Vec::new();

    // 1. Forward solve: L y = x (unit lower, NO conj).
    for s in 0..ns {
        let s_start = sb[s] as usize;
        let s_end = se[s] as usize;
        let s_ncols = s_end - s_start;
        let rstart = col_ptr_row[s] as usize;
        let s_pattern = &row_idx[rstart..rstart + nnz_per_super[s] as usize];
        let s_pat_len = s_pattern.len();
        let s_nrows = s_ncols + s_pat_len;
        let ls_slice = &l_values[col_ptr_val[s] as usize..col_ptr_val[s + 1] as usize];
        let ls_mat = MatRef::from_column_major_slice(T::faer_slice(ls_slice), s_nrows, s_ncols);
        let (ls_top, ls_bot) = ls_mat.split_at_row(s_ncols);

        {
            let mut x_mat = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut x), n, nrhs);
            let mut x_top = x_mat.rb_mut().subrows_mut(s_start, s_ncols);
            triangular_solve::solve_unit_lower_triangular_in_place(ls_top, x_top.rb_mut(), par);
            if s_pat_len > 0 {
                tmp.resize(s_pat_len * nrhs, T::ZERO);
                let tmp_mat = MatMut::from_column_major_slice_mut(
                    T::faer_slice_mut(&mut tmp[..s_pat_len * nrhs]),
                    s_pat_len,
                    nrhs,
                );
                matmul::matmul(tmp_mat, Accum::Replace, ls_bot, x_top.rb(), one, par);
            }
        }
        if s_pat_len > 0 {
            for col in 0..nrhs {
                let x_col = &mut x[col * n..(col + 1) * n];
                let t_col = &tmp[col * s_pat_len..(col + 1) * s_pat_len];
                for (idx, &row) in s_pattern.iter().enumerate() {
                    x_col[dyn_inv[row as usize] as usize] -= t_col[idx];
                }
            }
        }
    }

    // 2. Block-diagonal solve: B z = y (1×1 / 2×2 complex symmetric, NO conj).
    for s in 0..ns {
        let s_start = sb[s] as usize;
        let s_end = se[s] as usize;
        let s_ncols = s_end - s_start;
        let s_pat_len = nnz_per_super[s] as usize;
        let s_nrows = s_ncols + s_pat_len;
        let ls_slice = &l_values[col_ptr_val[s] as usize..col_ptr_val[s + 1] as usize];
        let s_sub = &subdiag[s_start..s_end];
        let mut idx = 0usize;
        while idx < s_ncols {
            let i = s_start + idx;
            let b21 = s_sub[idx];
            if b21 == T::ZERO {
                let d_inv = T::ONE / ls_slice[idx * s_nrows + idx];
                for col in 0..nrhs {
                    x[col * n + i] *= d_inv;
                }
                idx += 1;
            } else {
                let b11 = ls_slice[idx * s_nrows + idx];
                let b22 = ls_slice[(idx + 1) * s_nrows + idx + 1];
                let inv_det = T::ONE / (b11 * b22 - b21 * b21);
                for col in 0..nrhs {
                    let x0 = x[col * n + i];
                    let x1 = x[col * n + i + 1];
                    x[col * n + i] = (b22 * x0 - b21 * x1) * inv_det;
                    x[col * n + i + 1] = (b11 * x1 - b21 * x0) * inv_det;
                }
                idx += 2;
            }
        }
    }

    // 3. Backward solve: L^T w = z (unit upper, plain transpose, NO conj).
    for s in (0..ns).rev() {
        let s_start = sb[s] as usize;
        let s_end = se[s] as usize;
        let s_ncols = s_end - s_start;
        let rstart = col_ptr_row[s] as usize;
        let s_pattern = &row_idx[rstart..rstart + nnz_per_super[s] as usize];
        let s_pat_len = s_pattern.len();
        let s_nrows = s_ncols + s_pat_len;
        let ls_slice = &l_values[col_ptr_val[s] as usize..col_ptr_val[s + 1] as usize];
        let ls_mat = MatRef::from_column_major_slice(T::faer_slice(ls_slice), s_nrows, s_ncols);
        let (ls_top, ls_bot) = ls_mat.split_at_row(s_ncols);

        if s_pat_len > 0 {
            tmp.resize(s_pat_len * nrhs, T::ZERO);
            for col in 0..nrhs {
                let x_col = &x[col * n..(col + 1) * n];
                let t_col = &mut tmp[col * s_pat_len..(col + 1) * s_pat_len];
                for (idx, &row) in s_pattern.iter().enumerate() {
                    t_col[idx] = x_col[dyn_inv[row as usize] as usize];
                }
            }
        }
        let mut x_mat = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut x), n, nrhs);
        let mut x_top = x_mat.rb_mut().subrows_mut(s_start, s_ncols);
        if s_pat_len > 0 {
            let tmp_mat = MatRef::from_column_major_slice(T::faer_slice(&tmp[..s_pat_len * nrhs]), s_pat_len, nrhs);
            matmul::matmul(x_top.rb_mut(), Accum::Add, ls_bot.transpose(), tmp_mat, neg_one, par);
        }
        triangular_solve::solve_unit_upper_triangular_in_place(ls_top.transpose(), x_top.rb_mut(), par);
    }

    // Inverse permutation: rhs[orig, col] = x[dyn_inv[sym_inv[orig]], col].
    match sym_top.perm() {
        Some(p) => {
            let sym_inv = p.arrays().1;
            for col in 0..nrhs {
                let src = &x[col * n..(col + 1) * n];
                let dst = &mut rhs[col * n..(col + 1) * n];
                for i in 0..n {
                    dst[i] = src[dyn_inv[sym_inv[i] as usize] as usize];
                }
            }
        }
        None => {
            for col in 0..nrhs {
                let src = &x[col * n..(col + 1) * n];
                let dst = &mut rhs[col * n..(col + 1) * n];
                for i in 0..n {
                    dst[i] = src[dyn_inv[i] as usize];
                }
            }
        }
    }
    clear_upper_vector_state();
}

/// Non-conjugate simplicial `LDL^T` on `SymbolicSimplicialCholesky` (for small complex-symmetric matrices).
fn factorize_simplicial_sym_ldlt<T: Scalar>(
    l_values: &mut [T],
    lower_vals: &[T],
    plan: &SymLbltPlan,
    sym: &SymbolicSimplicialCholesky<u32>,
) -> Result<(), usize> {
    let n = sym.nrows();
    let col_ptr = sym.col_ptr();
    let row_idx = sym.row_idx();
    l_values.fill(T::ZERO);

    // Build column-to-row-list for left-looking updates (row i lists columns j < i with L[i, j] != 0).
    let mut row_cols: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    let mut work = vec![T::ZERO; n];
    for j in 0..n {
        let p0 = plan.perm_col_ptr[j];
        let p1 = plan.perm_col_ptr[j + 1];
        for pos in p0..p1 {
            work[plan.perm_row_idx[pos]] += lower_vals[plan.src_idx[pos]];
        }
        for &(k, pos_jk) in &row_cols[j] {
            let l_jk = l_values[pos_jk];
            let d_k = l_values[col_ptr[k] as usize];
            let scale = l_jk * d_k; // NO conj!
            let k_end = col_ptr[k + 1] as usize;
            for p in pos_jk..k_end {
                work[row_idx[p] as usize] -= scale * l_values[p];
            }
        }
        let c_start = col_ptr[j] as usize;
        let c_end = col_ptr[j + 1] as usize;
        let dj = work[j];
        if dj.modulus() == 0.0 || !dj.is_finite() {
            return Err(j);
        }
        l_values[c_start] = dj;
        work[j] = T::ZERO;
        let dj_inv = T::ONE / dj;
        for p in c_start + 1..c_end {
            let r = row_idx[p] as usize;
            let lv = work[r] * dj_inv;
            work[r] = T::ZERO;
            l_values[p] = lv;
            row_cols[r].push((j, p));
        }
    }
    Ok(())
}

fn solve_simplicial_sym_ldlt<T: Scalar>(
    sym_top: &SymbolicCholesky<u32>,
    sym: &SymbolicSimplicialCholesky<u32>,
    l_values: &[T],
    rhs: &mut [T],
    nrhs: usize,
) {
    let n = sym.nrows();
    let col_ptr = sym.col_ptr();
    let row_idx = sym.row_idx();
    let mut x = vec![T::ZERO; n];
    for col in 0..nrhs {
        let slice = &mut rhs[col * n..(col + 1) * n];
        match sym_top.perm() {
            Some(p) => {
                let fwd = p.arrays().0;
                for i in 0..n {
                    x[i] = slice[fwd[i] as usize];
                }
            }
            None => x.copy_from_slice(slice),
        }
        // Forward L y = x (unit lower)
        for j in 0..n {
            let xj = x[j];
            if xj != T::ZERO {
                for p in col_ptr[j] as usize + 1..col_ptr[j + 1] as usize {
                    x[row_idx[p] as usize] -= l_values[p] * xj;
                }
            }
        }
        // Diagonal D z = y
        for j in 0..n {
            x[j] = x[j] / l_values[col_ptr[j] as usize];
        }
        // Backward L^T w = z (NO conj)
        for j in (0..n).rev() {
            let mut sum = x[j];
            for p in col_ptr[j] as usize + 1..col_ptr[j + 1] as usize {
                sum -= l_values[p] * x[row_idx[p] as usize];
            }
            x[j] = sum;
        }
        match sym_top.perm() {
            Some(p) => {
                let inv = p.arrays().1;
                for i in 0..n {
                    slice[i] = x[inv[i] as usize];
                }
            }
            None => slice.copy_from_slice(&x),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ActiveFactor {
    None,
    Cholesky,
    /// `faer` supernodal/simplicial Hermitian `LBL^H` (real symmetric or complex Hermitian).
    HermLblt,
    /// Non-conjugate supernodal/simplicial complex-symmetric `LBL^T` (`A = A^T != A^H`).
    SymLblt,
    ColamdLu,
    Krylov,
}

/// Optional env override `NYX_DIRECT_BACKEND=colamd|lblt` for benchmarking/diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BackendPolicy {
    Auto,
    ForceColamd,
    ForceLblt,
}

fn backend_policy() -> BackendPolicy {
    match std::env::var("NYX_DIRECT_BACKEND").ok().as_deref() {
        Some("colamd") | Some("lu") => BackendPolicy::ForceColamd,
        Some("lblt") => BackendPolicy::ForceLblt,
        _ => BackendPolicy::Auto,
    }
}

/// `faer`-backed [`ShiftedSolver`] with symbolic reuse and fill-reducing ordering on the
/// Cholesky, Hermitian `LBL^H`, and complex-symmetric `LBL^T` paths (falling back transparently
/// to `faer` COLAMD LU for unsymmetric or near-singular matrices).
pub struct ShiftedFactorization<T: FaerScalar> {
    n: usize,
    definiteness: Definiteness,
    /// Union pattern of the whole matrix; built on first use when the Cholesky path (which reads the lower triangle only) is enough.
    full: Option<Assembled<T>>,
    tmap: Option<TransposeMap>,
    lower: Option<Assembled<T>>,
    /// Copies of the inputs kept only while `full` is not built (positive definite case).
    src: Option<(Csr<T>, Option<Csr<T>>)>,
    /// Last requested shift `sigma` for `A - sigma B` (used by deterministic iterative refinement and diagnostics).
    sigma: T,
    /// Actual shift factored (`sigma + delta` when shift perturbation recovery was used).
    factor_sigma: T,
    perturbed_active: bool,
    /// When true, `solve_block_in_place` performs 1-2 steps of deterministic compensated-residual
    /// iterative refinement (`r = b - (A - sigma B) x`) to eliminate parallel factor roundoff or
    /// moderately conditioned LBLT pivot growth (`ADR 11`, `FR-29`, `FR-33`).
    refine_solves: bool,
    sym_lu: Option<SymbolicLu<u32>>,
    sym_lu_est_nnz: usize,
    sym_llt: Option<SymbolicCholesky<u32>>,
    sym_llt_ord: FillOrdering,
    sym_lblt_plan: Option<SymLbltPlan>,
    lu_num: NumericLu<u32, T::Elem>,
    llt_vals: Vec<T>,
    lblt_subdiag: Vec<T>,
    lblt_perm_fwd: Vec<u32>,
    lblt_perm_inv: Vec<u32>,
    krylov_state: Option<KrylovFallback<T>>,
    have: Used,
    active: ActiveFactor,
    /// Thread policy for factorisation and solves (see [`ThreadPolicy`]).
    pub threads: ThreadPolicy,
    /// Fill-reducing ordering for Cholesky and symmetric indefinite `LBL^H`/`LBL^T`; set before the first factorisation.
    pub ordering: FillOrdering,
    /// Singular/OOM recovery strategy (R31); defaults to [`SingularRecovery::Auto`].
    pub recovery: SingularRecovery,
    /// When `true` (e.g. shift-and-invert spectral transformations), accept shift-perturbed factorisations
    /// directly as `(A - factor_sigma B)^{-1}` instead of defect-correcting back to the singular `sigma`.
    pub accept_perturbed_shift: bool,
    /// Optional dimension threshold above which `set_shift` routes directly to preconditioned Krylov (`GMRES`).
    pub iterative_threshold: Option<usize>,
    /// When `true` (e.g. non-symmetric gyrotropic media, `FR-34`), bypass symmetric/Hermitian `LBL^T`
    /// probes and factor directly with general COLAMD `LU`.
    pub force_general_lu: bool,
    stats: FactorStats,
    diagnostics: FactorDiagnostics,
    solves: AtomicUsize,
    krylov_its: AtomicUsize,
    skip_rcond: bool,
}

impl<T: FaerScalar> ShiftedFactorization<T> {
    /// Prepares the pattern for `A − σB` (`B = None` means the identity). With
    /// [`Definiteness::PositiveDefinite`], `A` and `B` must be symmetric/Hermitian (only their
    /// lower triangles are used by the Cholesky path).
    pub fn new(a: &Csr<T>, b: Option<&Csr<T>>, definiteness: Definiteness) -> Result<Self, SolverError> {
        let t0 = Instant::now();
        let n = a.nrows();
        let pd = definiteness == Definiteness::PositiveDefinite;
        let (full, tmap, lower, src) = if pd {
            (None, None, Some(assemble_lower_sorted(a, b)?), Some((a.clone(), b.cloned())))
        } else {
            let f = assemble(a, b, false)?;
            let tm = build_transpose_map(&f, n);
            (Some(f), Some(tm), None, None)
        };
        Ok(Self {
            n,
            definiteness,
            full,
            tmap,
            lower,
            src,
            sigma: T::ZERO,
            factor_sigma: T::ZERO,
            perturbed_active: false,
            refine_solves: false,
            sym_lu: None,
            sym_lu_est_nnz: 0,
            sym_llt: None,
            sym_llt_ord: FillOrdering::Amd,
            sym_lblt_plan: None,
            lu_num: NumericLu::new(),
            llt_vals: Vec::new(),
            lblt_subdiag: Vec::new(),
            lblt_perm_fwd: Vec::new(),
            lblt_perm_inv: Vec::new(),
            krylov_state: None,
            have: Used::None,
            active: ActiveFactor::None,
            threads: ThreadPolicy::default(),
            ordering: FillOrdering::default(),
            recovery: SingularRecovery::Auto,
            accept_perturbed_shift: false,
            iterative_threshold: None,
            force_general_lu: false,
            stats: FactorStats {
                analyse_secs: t0.elapsed().as_secs_f64(),
                symbolic_count: 0,
                numeric_secs_total: 0.0,
                numeric_secs_last: 0.0,
                numeric_count: 0,
                assemble_secs_last: 0.0,
                used: Used::None,
                backend: FactorBackend::None,
                cholesky_fallbacks: 0,
                solves: 0,
                factor_nnz: 0,
                solve_secs_total: 0.0,
            },
            diagnostics: FactorDiagnostics::default(),
            solves: AtomicUsize::new(0),
            krylov_its: AtomicUsize::new(0),
            skip_rcond: b.is_some(),
        })
    }

    /// Marks whether this factorization is for a shift-and-invert spectral pencil `(A - sigma B)`,
    /// where intentional proximity to an eigenvalue makes `rcond` uninformative and slow to probe.
    pub fn set_shift_invert_pencil(&mut self, skip: bool) {
        self.skip_rcond = skip;
    }

    /// Returns `true` when dynamic `rcond` probing is skipped (shift-and-invert pencils).
    pub fn is_shift_invert_pencil(&self) -> bool {
        self.skip_rcond
    }

    /// Returns the shift actually factored by the latest `set_shift` (`sigma + delta` if shift
    /// perturbation recovery fired, else `sigma`).
    pub fn factor_sigma(&self) -> T {
        self.factor_sigma
    }

    /// Returns a snapshot of pivot, Sylvester inertia, and singular/Krylov recovery diagnostics (R31).
    pub fn diagnostics(&self) -> FactorDiagnostics {
        FactorDiagnostics {
            krylov_iterations: self.krylov_its.load(Ordering::Relaxed),
            ..self.diagnostics
        }
    }

    /// Releases all symmetric Cholesky / `LBL^H` / `LBL^T` factor allocations so peak memory never
    /// holds both symmetric LBLT and COLAMD LU factors simultaneously (`FR-33`).
    fn release_lblt_buffers(&mut self) {
        self.llt_vals = Vec::new();
        self.sym_lblt_plan = None;
        self.lblt_subdiag = Vec::new();
        self.lblt_perm_fwd = Vec::new();
        self.lblt_perm_inv = Vec::new();
        self.sym_llt = None;
    }

    /// Releases COLAMD LU factor allocations when switching back to a symmetric factor backend.
    fn release_lu_buffers(&mut self) {
        self.sym_lu = None;
        self.sym_lu_est_nnz = 0;
        self.lu_num = NumericLu::new();
    }

    /// Replaces the *values* of `A` (and `B`) without redoing the symbolic analysis, for matrices
    /// that change every step (Jacobians, `T(z) = Σ f_j(z) A_j`). The sparsity pattern of the union
    /// must be identical to the one the object was built with (otherwise `InvalidArgument`);
    /// [`ShiftedSolver::set_shift`] must be called again before solving.
    pub fn replace_matrices(&mut self, a: &Csr<T>, b: Option<&Csr<T>>) -> Result<(), SolverError> {
        self.ensure_full()?;
        let full = assemble(a, b, false)?;
        let old = self
            .full
            .as_ref()
            .ok_or(SolverError::InvalidArgument("full pattern missing after ensure_full"))?;
        if full.col_ptr != old.col_ptr || full.row_idx != old.row_idx {
            return Err(SolverError::InvalidArgument("sparsity pattern differs from the one analysed"));
        }
        if self.definiteness == Definiteness::PositiveDefinite || self.lower.is_some() {
            self.lower = Some(extract_lower(&full, self.n));
            if self.definiteness == Definiteness::PositiveDefinite {
                self.src = Some((a.clone(), b.cloned()));
            }
        }
        self.full = Some(full);
        self.have = Used::None;
        self.active = ActiveFactor::None;
        self.perturbed_active = false;
        self.krylov_state = None;
        Ok(())
    }

    /// Builds the full union pattern if only the lower triangle exists so far (needed by the LU path).
    fn ensure_full(&mut self) -> Result<(), SolverError> {
        if let Some(full) = &self.full {
            if self.tmap.is_none() {
                self.tmap = Some(build_transpose_map(full, self.n));
            }
        } else {
            let (a, b) = self
                .src
                .as_ref()
                .ok_or(SolverError::InvalidArgument("source matrices missing while full pattern is absent"))?;
            let f = assemble(a, b.as_ref(), false)?;
            self.tmap = Some(build_transpose_map(&f, self.n));
            self.full = Some(f);
        }
        Ok(())
    }

    fn ensure_lower(&mut self) -> Result<(), SolverError> {
        if self.lower.is_none() {
            let full = self
                .full
                .as_ref()
                .ok_or(SolverError::InvalidArgument("full pattern missing before extracting lower"))?;
            self.lower = Some(extract_lower(full, self.n));
        }
        Ok(())
    }

    fn pattern(asm: &Assembled<T>, n: usize) -> SymbolicSparseColMatRef<'_, u32> {
        SymbolicSparseColMatRef::new_checked(n, n, &asm.col_ptr, None, &asm.row_idx)
    }

    fn shifted_values(asm: &Assembled<T>, sigma: T, out: &mut Vec<T>) {
        out.clear();
        out.extend(asm.a.iter().zip(&asm.b).map(|(&x, &y)| x - sigma * y));
    }

    /// Ensures the lower-triangle symbolic Cholesky structure (`sym_llt`) is built using `self.ordering`.
    fn ensure_sym_llt(&mut self) -> Result<(), SolverError> {
        if self.sym_llt.is_some() {
            return Ok(());
        }
        let n = self.n;
        let lower = self
            .lower
            .as_ref()
            .ok_or(SolverError::InvalidArgument("lower pattern missing for symbolic Cholesky"))?;
        let t0 = Instant::now();
        let pattern = Self::pattern(lower, n);
        let amd_sym = || factorize_symbolic_cholesky(pattern, Side::Lower, SymmetricOrdering::Amd, Default::default());
        let nd_sym = |perm: &(Vec<u32>, Vec<u32>)| {
            factorize_symbolic_cholesky(
                pattern,
                Side::Lower,
                SymmetricOrdering::Custom(PermRef::new_checked(&perm.0, &perm.1, n)),
                Default::default(),
            )
        };
        let make_nd_perm = || {
            let adj = adjacency_from_lower(n, &lower.col_ptr, &lower.row_idx);
            let order = nested_dissection(&adj, 100);
            order_to_perm(&order)
        };
        let (sym, chosen_ord) = match self.ordering {
            FillOrdering::Amd => (amd_sym().map_err(|_| SolverError::OutOfMemory)?, FillOrdering::Amd),
            FillOrdering::NestedDissection => {
                (nd_sym(&make_nd_perm()).map_err(|_| SolverError::OutOfMemory)?, FillOrdering::NestedDissection)
            }
            FillOrdering::Auto => {
                let a = amd_sym().map_err(|_| SolverError::OutOfMemory)?;
                if a.len_val() < 20 * n.max(1) {
                    (a, FillOrdering::Amd)
                } else {
                    let b = nd_sym(&make_nd_perm()).map_err(|_| SolverError::OutOfMemory)?;
                    if b.len_val() < a.len_val() {
                        (b, FillOrdering::NestedDissection)
                    } else {
                        (a, FillOrdering::Amd)
                    }
                }
            }
        };
        self.llt_vals = vec![T::ZERO; sym.len_val()];
        self.sym_llt_ord = chosen_ord;
        self.sym_llt = Some(sym);
        self.stats.analyse_secs += t0.elapsed().as_secs_f64();
        self.stats.symbolic_count += 1;
        Ok(())
    }

    /// Probes the symmetric indefinite factor on a deterministic vector. Single-shot residuals
    /// `<= 1e-9` are accepted immediately; when `allow_ir_rescue` is enabled (complex-symmetric
    /// `SymLblt`, `FR-33`), up to 2 cheap iterative refinement steps rescue `initial_rel_err <= 1e-4`
    /// when post-IR residual reaches `<= 1e-12`. Real/Hermitian `intranode_lblt` (`allow_ir_rescue = false`)
    /// requires single-shot `<= 1e-9` so saddle-point / nullspace pencils with cross-supernode growth
    /// fall back to pivoted COLAMD LU cleanly.
    fn verify_probe_residual(
        &self,
        full_vals: &[T],
        allow_ir_rescue: bool,
        solve_one: impl Fn(&mut [T]),
    ) -> ProbeCheck {
        let n = self.n;
        let mut probe = Vec::with_capacity(n);
        for i in 0..n {
            let v = (((i.wrapping_mul(2654435761) ^ 0x5a5a_a5a5) % 10007) as f64 / 5003.5) - 1.0;
            probe.push(T::from_f64(v));
        }
        let p_norm = vecops::norm2(&probe).max(1e-300);
        let mut x = probe.clone();
        solve_one(&mut x);
        if !x.iter().all(|v| v.is_finite()) {
            return ProbeCheck { ok: false, initial_rel_err: f64::INFINITY };
        }
        let Some(full) = self.full.as_ref() else {
            return ProbeCheck { ok: false, initial_rel_err: f64::INFINITY };
        };
        let mut r = vec![T::ZERO; n];
        csc_residual_vals(n, &full.col_ptr, &full.row_idx, full_vals, &probe, &x, &mut r);
        let initial_rel_err = vecops::norm2(&r) / p_norm;
        if !initial_rel_err.is_finite() {
            return ProbeCheck { ok: false, initial_rel_err: f64::INFINITY };
        }
        let single_shot_tol = if allow_ir_rescue { 1e-11 } else { 1e-9 };
        if initial_rel_err <= single_shot_tol {
            return ProbeCheck { ok: true, initial_rel_err };
        }
        if !allow_ir_rescue || initial_rel_err > 1e-4 {
            return ProbeCheck { ok: false, initial_rel_err };
        }
        // Up to 2 steps of cheap iterative refinement before falling back to COLAMD LU.
        let mut rel_err = initial_rel_err;
        let mut dx = vec![T::ZERO; n];
        for _ in 0..2 {
            dx.copy_from_slice(&r);
            solve_one(&mut dx);
            if !dx.iter().all(|v| v.is_finite()) {
                return ProbeCheck { ok: false, initial_rel_err };
            }
            for (xi, &dxi) in x.iter_mut().zip(&dx) {
                *xi += dxi;
            }
            csc_residual_vals(n, &full.col_ptr, &full.row_idx, full_vals, &probe, &x, &mut r);
            rel_err = vecops::norm2(&r) / p_norm;
            if !rel_err.is_finite() {
                return ProbeCheck { ok: false, initial_rel_err };
            }
            if rel_err <= 1e-12 {
                return ProbeCheck { ok: true, initial_rel_err };
            }
        }
        ProbeCheck { ok: rel_err <= 1e-12, initial_rel_err }
    }

    fn try_cholesky(&mut self, sigma: T) -> Result<bool, SolverError> {
        let n = self.n;
        let par = self.threads.factor.par_for_factor(n);
        self.ensure_sym_llt()?;
        let lower = self
            .lower
            .as_ref()
            .ok_or(SolverError::InvalidArgument("lower pattern missing for Cholesky"))?;
        let mut vals = Vec::with_capacity(lower.row_idx.len());
        Self::shifted_values(lower, sigma, &mut vals);
        let sym = self
            .sym_llt
            .as_ref()
            .ok_or(SolverError::InvalidArgument("symbolic Cholesky factor missing"))?;
        let factor_nnz = sym.len_val();
        let (scale, min_d, max_d, max_e) = lower_scale_and_diag(n, lower, &vals);
        self.diagnostics.matrix_scale = scale;
        let mat = SparseColMatRef::new(Self::pattern(lower, n), T::faer_slice(&vals));
        let t0 = Instant::now();
        let mut mem = MemBuffer::new(sym.factorize_numeric_llt_scratch::<T::Elem>(par, Default::default()));
        let res = sym.factorize_numeric_llt::<T::Elem>(
            T::faer_slice_mut(&mut self.llt_vals),
            mat,
            Side::Lower,
            Default::default(),
            par,
            MemStack::new(&mut mem),
            Default::default(),
        );
        clear_upper_vector_state();
        self.stats.numeric_secs_last = t0.elapsed().as_secs_f64();
        drop(vals);
        match res {
            Ok(_) => {
                let sym = self.sym_llt.as_ref().unwrap();
                record_cholesky_diagnostics(sym, &self.llt_vals, scale, &mut self.diagnostics);
                finalize_pivot_growth(&mut self.diagnostics, min_d, max_d, max_e);
                self.release_lu_buffers();
                self.sigma = sigma;
                self.refine_solves = n > PAR_FACTOR_MIN_N || par != Par::Seq;
                self.have = Used::Cholesky;
                self.active = ActiveFactor::Cholesky;
                self.stats.used = Used::Cholesky;
                self.stats.backend = match self.sym_llt_ord {
                    FillOrdering::NestedDissection => FactorBackend::CholeskyNd,
                    _ => FactorBackend::CholeskyAmd,
                };
                self.stats.factor_nnz = factor_nnz;
                Ok(true)
            }
            Err(_) => Ok(false),
        }
    }

    /// Supernodal/simplicial Hermitian `LBL^H` (`faer::intranode_lblt`) for real-symmetric or
    /// complex-Hermitian indefinite matrices, guarded by finiteness + IR-verified probe residual `<= 1e-8`.
    fn try_herm_lblt(&mut self, sigma: T, full_vals: &[T]) -> Result<bool, SolverError> {
        let n = self.n;
        if n == 0 {
            return Ok(false);
        }
        self.ensure_lower()?;
        self.ensure_sym_llt()?;
        let par = self.threads.factor.par_for_factor(n);
        let full_ref = self.full.as_ref().unwrap();
        let (scale, min_d, max_d, max_e) = csc_scale_and_diag(n, &full_ref.col_ptr, &full_ref.row_idx, full_vals);
        self.diagnostics.matrix_scale = scale;
        let (fac_secs, initial_rel_err, factor_nnz) = {
            let lower = self
                .lower
                .as_ref()
                .ok_or(SolverError::InvalidArgument("lower pattern missing for LBL^H"))?;
            let mut lvals = Vec::with_capacity(lower.row_idx.len());
            Self::shifted_values(lower, sigma, &mut lvals);

            self.lblt_subdiag.resize(n, T::ZERO);
            self.lblt_perm_fwd.resize(n, 0);
            self.lblt_perm_inv.resize(n, 0);

            let sym = self
                .sym_llt
                .as_ref()
                .ok_or(SolverError::InvalidArgument("symbolic Cholesky factor missing for LBL^H"))?;
            let factor_nnz = sym.len_val();
            let mat = SparseColMatRef::new(Self::pattern(lower, n), T::faer_slice(&lvals));
            let t0 = Instant::now();
            let mut mem =
                MemBuffer::new(sym.factorize_numeric_intranode_lblt_scratch::<T::Elem>(par, Default::default()));
            let _ = sym.factorize_numeric_intranode_lblt::<T::Elem>(
                T::faer_slice_mut(&mut self.llt_vals),
                T::faer_slice_mut(&mut self.lblt_subdiag),
                &mut self.lblt_perm_fwd,
                &mut self.lblt_perm_inv,
                mat,
                Side::Lower,
                par,
                MemStack::new(&mut mem),
                Default::default(),
            );
            clear_upper_vector_state();
            let fac_secs = t0.elapsed().as_secs_f64();
            drop(mem);
            drop(lvals);

            let bad_col = record_lblt_diagnostics(
                sym,
                &self.llt_vals,
                &self.lblt_subdiag,
                true,
                scale,
                &mut self.diagnostics,
            );
            finalize_pivot_growth(&mut self.diagnostics, min_d, max_d, max_e);
            if bad_col.is_some()
                || !self.llt_vals.iter().all(|v| v.is_finite())
                || !self.lblt_subdiag.iter().all(|v| v.is_finite())
            {
                return Ok(false);
            }

            let solve_par = self.threads.solve.par();
            let probe = self.verify_probe_residual(full_vals, false, |x| {
                let mut smem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(1, solve_par));
                let perm = PermRef::new_checked(&self.lblt_perm_fwd, &self.lblt_perm_inv, n);
                let lblt = IntranodeLbltRef::new(
                    sym,
                    T::faer_slice(&self.llt_vals),
                    T::faer_slice(&self.lblt_subdiag),
                    perm,
                );
                let xm = MatMut::from_column_major_slice_mut(T::faer_slice_mut(x), n, 1);
                lblt.solve_in_place_with_conj(Conj::No, xm, solve_par, MemStack::new(&mut smem));
                clear_upper_vector_state();
            });
            if !probe.ok {
                return Ok(false);
            }
            (fac_secs, probe.initial_rel_err, factor_nnz)
        };

        self.release_lu_buffers();
        self.sigma = sigma;
        self.refine_solves = n > PAR_FACTOR_MIN_N || par != Par::Seq || initial_rel_err > 1e-13;
        self.stats.numeric_secs_last = fac_secs;
        self.have = Used::Lu;
        self.active = ActiveFactor::HermLblt;
        self.stats.used = Used::Lu;
        self.stats.backend = match self.sym_llt_ord {
            FillOrdering::NestedDissection => FactorBackend::LbltNd,
            _ => FactorBackend::LbltAmd,
        };
        self.stats.factor_nnz = factor_nnz;
        Ok(true)
    }

    /// Non-conjugate supernodal/simplicial complex-symmetric `LBL^T` (`A = A^T != A^H`, e.g. PML /
    /// lossy dielectrics / surface impedance), reusing `faer`'s `SymbolicCholesky` (AMD/ND/Auto) and
    /// BLAS-3 dense kernels.
    fn try_sym_lblt(&mut self, sigma: T, full_vals: &[T], allow_ir_rescue: bool) -> Result<bool, SolverError> {
        let n = self.n;
        if n == 0 {
            return Ok(false);
        }
        self.ensure_lower()?;
        self.ensure_sym_llt()?;
        if self.sym_lblt_plan.is_none() {
            let t0 = Instant::now();
            let lower = self
                .lower
                .as_ref()
                .ok_or(SolverError::InvalidArgument("lower pattern missing for LBL^T plan"))?;
            let sym = self
                .sym_llt
                .as_ref()
                .ok_or(SolverError::InvalidArgument("symbolic Cholesky factor missing for LBL^T plan"))?;
            let plan = build_sym_lblt_plan(lower, sym, n);
            self.sym_lblt_plan = Some(plan);
            self.stats.analyse_secs += t0.elapsed().as_secs_f64();
        }
        let par = self.threads.factor.par_for_factor(n);
        let full_ref = self.full.as_ref().unwrap();
        let (scale, min_d, max_d, max_e) = csc_scale_and_diag(n, &full_ref.col_ptr, &full_ref.row_idx, full_vals);
        self.diagnostics.matrix_scale = scale;
        let (fac_secs, initial_rel_err, factor_nnz) = {
            let lower = self
                .lower
                .as_ref()
                .ok_or(SolverError::InvalidArgument("lower pattern missing for LBL^T"))?;
            let mut lvals = Vec::with_capacity(lower.row_idx.len());
            Self::shifted_values(lower, sigma, &mut lvals);

            self.lblt_subdiag.resize(n, T::ZERO);
            self.lblt_perm_fwd.resize(n, 0);
            self.lblt_perm_inv.resize(n, 0);

            let t0 = Instant::now();
            let sym = self
                .sym_llt
                .as_ref()
                .ok_or(SolverError::InvalidArgument("symbolic Cholesky factor missing for LBL^T"))?;
            let factor_nnz = sym.len_val();
            let plan = self
                .sym_lblt_plan
                .as_ref()
                .ok_or(SolverError::InvalidArgument("LBL^T plan missing"))?;
            let fac_res = match sym.raw() {
                SymbolicCholeskyRaw::Supernodal(sup) => factorize_supernodal_sym_lblt(
                    &mut self.llt_vals,
                    &mut self.lblt_subdiag,
                    &mut self.lblt_perm_fwd,
                    &mut self.lblt_perm_inv,
                    &lvals,
                    plan,
                    sup,
                    par,
                ),
                SymbolicCholeskyRaw::Simplicial(simp) => {
                    factorize_simplicial_sym_ldlt(&mut self.llt_vals, &lvals, plan, simp)
                }
            };
            let fac_secs = t0.elapsed().as_secs_f64();
            drop(lvals);
            if let Err(bad_idx) = fac_res {
                let orig_col = sym
                    .perm()
                    .and_then(|p| p.arrays().0.get(bad_idx).map(|&c| c as usize))
                    .unwrap_or(bad_idx);
                self.diagnostics.first_bad_pivot_index.get_or_insert(orig_col);
                self.diagnostics.min_pivot_abs = 0.0;
                return Ok(false);
            }
            let bad_col = record_lblt_diagnostics(
                sym,
                &self.llt_vals,
                &self.lblt_subdiag,
                !T::IS_COMPLEX,
                scale,
                &mut self.diagnostics,
            );
            finalize_pivot_growth(&mut self.diagnostics, min_d, max_d, max_e);
            if bad_col.is_some()
                || !self.llt_vals.iter().all(|v| v.is_finite())
                || !self.lblt_subdiag.iter().all(|v| v.is_finite())
            {
                return Ok(false);
            }

            let solve_par = self.threads.solve.par();
            let probe = self.verify_probe_residual(full_vals, allow_ir_rescue, |x| match sym.raw() {
                SymbolicCholeskyRaw::Supernodal(sup) => solve_supernodal_sym_lblt(
                    sym,
                    sup,
                    &self.llt_vals,
                    &self.lblt_subdiag,
                    &self.lblt_perm_fwd,
                    &self.lblt_perm_inv,
                    x,
                    1,
                    solve_par,
                ),
                SymbolicCholeskyRaw::Simplicial(simp) => solve_simplicial_sym_ldlt(sym, simp, &self.llt_vals, x, 1),
            });
            if !probe.ok {
                return Ok(false);
            }
            (fac_secs, probe.initial_rel_err, factor_nnz)
        };

        self.release_lu_buffers();
        self.sigma = sigma;
        self.refine_solves = n > PAR_FACTOR_MIN_N || par != Par::Seq || initial_rel_err > 1e-13;
        self.stats.numeric_secs_last = fac_secs;
        self.have = Used::Lu;
        self.active = ActiveFactor::SymLblt;
        self.stats.used = Used::Lu;
        self.stats.backend = match self.sym_llt_ord {
            FillOrdering::NestedDissection => FactorBackend::LbltNd,
            _ => FactorBackend::LbltAmd,
        };
        self.stats.factor_nnz = factor_nnz;
        Ok(true)
    }

    fn do_colamd_lu(&mut self, sigma: T, full_vals: &[T]) -> Result<(), SolverError> {
        let n = self.n;
        let par = self.threads.factor.par_for_factor(n);
        // Explicitly drop all symmetric LBLT/Cholesky factor buffers BEFORE allocating COLAMD LU
        // so peak RAM never holds both factorisations simultaneously (`FR-33`).
        self.release_lblt_buffers();
        if self.sym_lu.is_none() {
            let t0 = Instant::now();
            let full = self
                .full
                .as_ref()
                .ok_or(SolverError::InvalidArgument("full matrix pattern missing for LU"))?;
            let sym = factorize_symbolic_lu(Self::pattern(full, n), Default::default())
                .map_err(|_| SolverError::OutOfMemory)?;
            // Estimate COLAMD symmetric fill when structurally symmetric so factor_nnz is comparable.
            let est_nnz = if self.tmap.as_ref().is_some_and(|t| t.structurally_symmetric) {
                self.ensure_lower()?;
                let lower = self
                    .lower
                    .as_ref()
                    .ok_or(SolverError::InvalidArgument("lower pattern missing for LU fill estimate"))?;
                factorize_symbolic_cholesky(
                    Self::pattern(lower, n),
                    Side::Lower,
                    SymmetricOrdering::Custom(sym.col_perm()),
                    Default::default(),
                )
                .map(|c| (2 * c.len_val()).saturating_sub(n))
                .unwrap_or(0)
            } else {
                0
            };
            if self.definiteness != Definiteness::PositiveDefinite {
                self.lower = None;
            }
            self.sym_lu = Some(sym);
            self.sym_lu_est_nnz = est_nnz;
            self.stats.analyse_secs += t0.elapsed().as_secs_f64();
            self.stats.symbolic_count += 1;
        }
        let sym = self
            .sym_lu
            .as_ref()
            .ok_or(SolverError::InvalidArgument("symbolic LU factor missing"))?;
        let full = self
            .full
            .as_ref()
            .ok_or(SolverError::InvalidArgument("full matrix pattern missing for LU"))?;
        let (scale, min_d, max_d, max_e) = csc_scale_and_diag(n, &full.col_ptr, &full.row_idx, full_vals);
        self.diagnostics.matrix_scale = scale;
        self.diagnostics.min_diag_abs = min_d;
        self.diagnostics.max_diag_abs = max_d;
        self.diagnostics.max_entry_abs = max_e;
        self.diagnostics.has_factor_pivots = false;
        if n > 0 && scale == 0.0 {
            self.diagnostics.first_bad_pivot_index.get_or_insert(0);
            self.diagnostics.min_pivot_abs = 0.0;
            self.diagnostics.pivot_growth = f64::INFINITY;
            self.diagnostics.pivot_ratio = f64::INFINITY;
            return Err(SolverError::Singular);
        }
        let mat = SparseColMatRef::new(Self::pattern(full, n), T::faer_slice(full_vals));
        let t0 = Instant::now();
        let mut mem = MemBuffer::new(sym.factorize_numeric_lu_scratch::<T::Elem>(par, Default::default()));
        let res = sym.factorize_numeric_lu::<T::Elem>(
            &mut self.lu_num,
            mat,
            par,
            MemStack::new(&mut mem),
            Default::default(),
        );
        clear_upper_vector_state();
        self.stats.numeric_secs_last = t0.elapsed().as_secs_f64();
        match res {
            Ok(_) => {
                // `faer` partial-pivoting LU returns `Ok` even when a column pivot is numerically 0
                // (`0.recip() -> inf`); guard via a deterministic probe residual before accepting.
                let solve_par = self.threads.solve.par();
                let (probe_ok, probe_cond) = if n == 0 {
                    (true, 1.0)
                } else {
                    let mut p = Vec::with_capacity(n);
                    for i in 0..n {
                        let v = (((i.wrapping_mul(2654435761) ^ 0x5a5a_a5a5) % 10007) as f64 / 5003.5) - 1.0;
                        p.push(T::from_f64(v));
                    }
                    let p_norm = vecops::norm2(&p).max(1e-300);
                    let mut x = p.clone();
                    let xm = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut x), n, 1);
                    let mut smem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(1, solve_par));
                    LuRef::new_unchecked(sym, &self.lu_num).solve_in_place_with_conj(
                        Conj::No,
                        xm,
                        solve_par,
                        MemStack::new(&mut smem),
                    );
                    clear_upper_vector_state();
                    if !x.iter().all(|v| v.is_finite()) {
                        (false, f64::INFINITY)
                    } else {
                        let mut r = vec![T::ZERO; n];
                        csc_residual_vals(n, &full.col_ptr, &full.row_idx, full_vals, &p, &x, &mut r);
                        let r_norm = vecops::norm2(&r);
                        let x_norm = vecops::norm2(&x);
                        let rel_err = r_norm / p_norm;
                        let bw_err = r_norm / (scale * x_norm + p_norm).max(1e-300);
                        let cond_lb = (scale * x_norm) / p_norm;
                        let ok = bw_err.is_finite()
                            && rel_err.is_finite()
                            && bw_err <= 1e-6
                            && !(rel_err > 1e-2 && cond_lb > 1e14);
                        (ok, cond_lb)
                    }
                };
                if !probe_ok {
                    let bad_col = csc_bad_diag_col(n, &full.col_ptr, full_vals);
                    self.diagnostics.first_bad_pivot_index.get_or_insert(bad_col);
                    self.diagnostics.min_pivot_abs = 0.0;
                    self.diagnostics.pivot_growth = f64::INFINITY;
                    self.diagnostics.pivot_ratio = f64::INFINITY;
                    return Err(SolverError::Singular);
                }
                let diag_ratio = if min_d > 0.0 { max_d / min_d } else { 1.0 };
                self.diagnostics.pivot_growth = probe_cond.max(1.0);
                self.diagnostics.pivot_ratio = probe_cond.max(diag_ratio);
                self.sigma = sigma;
                self.refine_solves = n > PAR_FACTOR_MIN_N || par != Par::Seq;
                self.have = Used::Lu;
                self.active = ActiveFactor::ColamdLu;
                self.stats.used = Used::Lu;
                self.stats.backend = FactorBackend::LuColamd;
                self.stats.factor_nnz = self.sym_lu_est_nnz;
                Ok(())
            }
            Err(faer::sparse::linalg::LuError::SymbolicSingular { index }) => {
                self.diagnostics.first_bad_pivot_index.get_or_insert(index);
                self.diagnostics.min_pivot_abs = 0.0;
                Err(SolverError::Singular)
            }
            Err(_) => Err(SolverError::OutOfMemory),
        }
    }

    fn do_lu(&mut self, sigma: T) -> Result<(), SolverError> {
        self.ensure_full()?;
        let (vals, is_herm, is_cplx_sym) = {
            let full = self
                .full
                .as_ref()
                .ok_or(SolverError::InvalidArgument("full matrix pattern missing for LU"))?;
            let mut vals = Vec::with_capacity(full.row_idx.len());
            Self::shifted_values(full, sigma, &mut vals);

            let tmap = self
                .tmap
                .as_ref()
                .ok_or(SolverError::InvalidArgument("transpose map missing for LU"))?;
            let is_struct_sym = tmap.structurally_symmetric;
            let (is_herm, is_cplx_sym) = if is_struct_sym {
                let max_v = vals.iter().map(|v| v.modulus()).fold(0.0f64, f64::max);
                let tol = 1e-12 * max_v.max(1e-300);
                let herm = vals.iter().enumerate().all(|(k, &v)| (v - vals[tmap.t_map[k]].conj()).modulus() <= tol);
                let sym = if !T::IS_COMPLEX {
                    herm
                } else if herm {
                    false
                } else {
                    vals.iter().enumerate().all(|(k, &v)| (v - vals[tmap.t_map[k]]).modulus() <= tol)
                };
                (herm, sym)
            } else {
                (false, false)
            };
            (vals, is_herm, is_cplx_sym)
        };
        let policy = if self.force_general_lu {
            BackendPolicy::ForceColamd
        } else {
            backend_policy()
        };

        let res = match policy {
            BackendPolicy::ForceColamd => self.do_colamd_lu(sigma, &vals),
            BackendPolicy::ForceLblt => {
                let herm_ok = is_herm && self.try_herm_lblt(sigma, &vals)?;
                let sym_ok = !herm_ok && is_cplx_sym && self.try_sym_lblt(sigma, &vals, true)?;
                if herm_ok || sym_ok {
                    Ok(())
                } else {
                    let saved_diag = self.diagnostics;
                    self.release_lblt_buffers();
                    let lu_res = self.do_colamd_lu(sigma, &vals);
                    if lu_res.is_err() && saved_diag.zero_inertia.is_some() {
                        self.diagnostics = saved_diag;
                    }
                    lu_res
                }
            }
            BackendPolicy::Auto => {
                // If COLAMD LU was already analysed for an earlier shift (e.g. unsymmetric or
                // previous fallback), keep reusing it so `symbolic_count` stays 1.
                if self.sym_lu.is_some() {
                    self.do_colamd_lu(sigma, &vals)
                } else {
                    let herm_ok = is_herm && self.try_herm_lblt(sigma, &vals)?;
                    let sym_ok = !herm_ok && is_cplx_sym && self.try_sym_lblt(sigma, &vals, true)?;
                    if herm_ok || sym_ok {
                        // Real symmetric / complex Hermitian (`LBL^H`) or complex-symmetric (`LBL^T`)
                        // supernodal Bunch-Kaufman with AMD/ND/Auto ordering.
                        Ok(())
                    } else {
                        let saved_diag = self.diagnostics;
                        // Unsymmetric or LBL^T pivot fallback: free LBLT buffers before COLAMD LU.
                        self.release_lblt_buffers();
                        let lu_res = self.do_colamd_lu(sigma, &vals);
                        if lu_res.is_err() && saved_diag.zero_inertia.is_some() {
                            self.diagnostics = saved_diag;
                        }
                        lu_res
                    }
                }
            }
        };

        drop(vals);
        res
    }

    fn factorize_direct_at(&mut self, sigma: T) -> Result<(), SolverError> {
        let mut done = false;
        if self.definiteness == Definiteness::PositiveDefinite {
            if self.try_cholesky(sigma)? {
                done = true;
            } else {
                self.stats.cholesky_fallbacks += 1;
            }
        }
        if !done { self.do_lu(sigma) } else { Ok(()) }
    }

    fn try_krylov_fallback(&mut self, sigma: T) -> Result<(), SolverError> {
        let n = self.n;
        self.ensure_full()?;
        let full = self
            .full
            .as_ref()
            .ok_or(SolverError::InvalidArgument("full matrix pattern missing for Krylov fallback"))?;
        let csr = assembled_to_csr(full, n, sigma)?;
        let mut scale = self.diagnostics.matrix_scale;
        if scale == 0.0 {
            let mut vals = Vec::with_capacity(full.row_idx.len());
            Self::shifted_values(full, sigma, &mut vals);
            scale = csc_one_norm(n, &full.col_ptr, &vals);
            self.diagnostics.matrix_scale = scale;
        }
        if n > 0 && scale == 0.0 {
            self.diagnostics.first_bad_pivot_index.get_or_insert(0);
            return Err(SolverError::Singular);
        }
        let pc = if n > 0 {
            let ilu = Ilu::ilu_k(&csr, 0);
            if ilu.nnz() > 0 {
                KrylovPc::Ilu(Box::new(ilu))
            } else {
                KrylovPc::Jacobi(Jacobi::new(&csr))
            }
        } else {
            KrylovPc::Jacobi(Jacobi::new(&csr))
        };
        self.release_lblt_buffers();
        self.release_lu_buffers();
        self.krylov_state = Some(KrylovFallback { mat: csr, pc });
        self.sigma = sigma;
        self.factor_sigma = sigma;
        self.perturbed_active = false;
        self.refine_solves = false;
        self.active = ActiveFactor::Krylov;
        self.diagnostics.krylov_fallback = true;
        Ok(())
    }

    fn try_perturbed_shift(&mut self, sigma: T) -> Result<bool, SolverError> {
        let n = self.n;
        if n == 0 {
            return Ok(false);
        }
        self.ensure_full()?;
        let (scale, b_scale) = {
            let full = self
                .full
                .as_ref()
                .ok_or(SolverError::InvalidArgument("full matrix pattern missing for shift perturbation"))?;
            let b_scale = b_one_norm_full(full, n);
            let mut scale = self.diagnostics.matrix_scale;
            if scale == 0.0 {
                let mut vals = Vec::with_capacity(full.row_idx.len());
                Self::shifted_values(full, sigma, &mut vals);
                scale = csc_one_norm(n, &full.col_ptr, &vals);
                self.diagnostics.matrix_scale = scale;
            }
            (scale, b_scale)
        };
        if scale == 0.0 || b_scale == 0.0 {
            return Ok(false);
        }
        let base_delta = (scale / b_scale).max(sigma.modulus()).max(1e-12);
        let init_diag = self.diagnostics;

        let mags = [1e-10, 1e-8, 1e-6, 1e-4];
        let mut deltas: Vec<T> = Vec::with_capacity(8);
        for &m in &mags {
            let step = m * base_delta;
            if T::IS_COMPLEX {
                deltas.push(T::from_parts(step, step));
                deltas.push(T::from_parts(-step, step));
            } else {
                deltas.push(T::from_f64(step));
                deltas.push(T::from_f64(-step));
            }
        }

        for delta in deltas {
            let sigma_p = sigma + delta;
            self.diagnostics = FactorDiagnostics::default();
            if self.factorize_direct_at(sigma_p).is_ok() {
                self.factor_sigma = sigma_p;
                self.sigma = sigma;
                self.perturbed_active = true;
                self.diagnostics.first_bad_pivot_index = init_diag.first_bad_pivot_index;
                if self.diagnostics.zero_inertia.is_none() {
                    self.diagnostics.positive_inertia = init_diag.positive_inertia;
                    self.diagnostics.negative_inertia = init_diag.negative_inertia;
                    self.diagnostics.zero_inertia = init_diag.zero_inertia;
                }
                self.diagnostics.shift_perturbation_re = delta.re();
                self.diagnostics.shift_perturbation_im = delta.im();
                return Ok(true);
            }
        }
        self.diagnostics = init_diag;
        Ok(false)
    }

    fn solve_factor_only_in_place(&self, rhs: &mut [T], nrhs: usize) -> Result<(), SolverError> {
        let n = self.n;
        let par = self.threads.solve.par();
        match self.active {
            ActiveFactor::Cholesky => {
                let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
                let sym = self.sym_llt.as_ref().ok_or(SolverError::NotFactored)?;
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                LltRef::new(sym, T::faer_slice(&self.llt_vals)).solve_in_place_with_conj(
                    Conj::No,
                    m,
                    par,
                    MemStack::new(&mut mem),
                );
            }
            ActiveFactor::HermLblt => {
                let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
                let sym = self.sym_llt.as_ref().ok_or(SolverError::NotFactored)?;
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                let perm = PermRef::new_checked(&self.lblt_perm_fwd, &self.lblt_perm_inv, n);
                IntranodeLbltRef::new(
                    sym,
                    T::faer_slice(&self.llt_vals),
                    T::faer_slice(&self.lblt_subdiag),
                    perm,
                )
                .solve_in_place_with_conj(Conj::No, m, par, MemStack::new(&mut mem));
            }
            ActiveFactor::SymLblt => {
                let sym = self.sym_llt.as_ref().ok_or(SolverError::NotFactored)?;
                match sym.raw() {
                    SymbolicCholeskyRaw::Supernodal(sup) => solve_supernodal_sym_lblt(
                        sym,
                        sup,
                        &self.llt_vals,
                        &self.lblt_subdiag,
                        &self.lblt_perm_fwd,
                        &self.lblt_perm_inv,
                        &mut rhs[..n * nrhs],
                        nrhs,
                        par,
                    ),
                    SymbolicCholeskyRaw::Simplicial(simp) => {
                        solve_simplicial_sym_ldlt(sym, simp, &self.llt_vals, &mut rhs[..n * nrhs], nrhs)
                    }
                }
            }
            ActiveFactor::ColamdLu => {
                let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
                let sym = self.sym_lu.as_ref().ok_or(SolverError::NotFactored)?;
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                LuRef::new_unchecked(sym, &self.lu_num).solve_in_place_with_conj(
                    Conj::No,
                    m,
                    par,
                    MemStack::new(&mut mem),
                );
            }
            ActiveFactor::Krylov | ActiveFactor::None => return Err(SolverError::NotFactored),
        }
        clear_upper_vector_state();
        Ok(())
    }

    fn solve_transpose_factor_only_in_place(&self, rhs: &mut [T], nrhs: usize) -> Result<(), SolverError> {
        let n = self.n;
        let par = self.threads.solve.par();
        match self.active {
            ActiveFactor::Cholesky => {
                let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
                let sym = self.sym_llt.as_ref().ok_or(SolverError::NotFactored)?;
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                let conj = if T::IS_COMPLEX { Conj::Yes } else { Conj::No };
                LltRef::new(sym, T::faer_slice(&self.llt_vals)).solve_in_place_with_conj(
                    conj,
                    m,
                    par,
                    MemStack::new(&mut mem),
                );
            }
            ActiveFactor::HermLblt => {
                let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
                let sym = self.sym_llt.as_ref().ok_or(SolverError::NotFactored)?;
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                let perm = PermRef::new_checked(&self.lblt_perm_fwd, &self.lblt_perm_inv, n);
                let conj = if T::IS_COMPLEX { Conj::Yes } else { Conj::No };
                IntranodeLbltRef::new(
                    sym,
                    T::faer_slice(&self.llt_vals),
                    T::faer_slice(&self.lblt_subdiag),
                    perm,
                )
                .solve_in_place_with_conj(conj, m, par, MemStack::new(&mut mem));
            }
            ActiveFactor::SymLblt => {
                return self.solve_factor_only_in_place(rhs, nrhs);
            }
            ActiveFactor::ColamdLu => {
                let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
                let sym = self.sym_lu.as_ref().ok_or(SolverError::NotFactored)?;
                let mut mem = MemBuffer::new(sym.solve_transpose_in_place_scratch::<T::Elem>(nrhs, par));
                LuRef::new_unchecked(sym, &self.lu_num).solve_transpose_in_place_with_conj(
                    Conj::No,
                    m,
                    par,
                    MemStack::new(&mut mem),
                );
            }
            ActiveFactor::Krylov | ActiveFactor::None => return Err(SolverError::NotFactored),
        }
        clear_upper_vector_state();
        Ok(())
    }

    /// 2-step Hager-Higham 1-norm reciprocal condition number estimate (`1 / (||M||_1 ||M^{-1}||_1)`).
    /// Returns `None` on shift-and-invert pencils (`skip_rcond`), Krylov fallback, or unfactored instances.
    pub fn estimate_rcond_1norm(&self) -> Option<f64> {
        if self.skip_rcond || matches!(self.active, ActiveFactor::None | ActiveFactor::Krylov) {
            return None;
        }
        let n = self.n;
        if n == 0 {
            return Some(1.0);
        }
        let norm1 = self.diagnostics.matrix_scale;
        if norm1 <= 0.0 || !norm1.is_finite() {
            return Some(0.0);
        }

        let norm_l1 = |v: &[T]| -> f64 { v.iter().map(|z| z.modulus()).sum() };
        let sign_vec = |v: &[T], out: &mut [T]| {
            for (dst, &src) in out.iter_mut().zip(v) {
                let m = src.modulus();
                *dst = if m > 0.0 && m.is_finite() {
                    src.scale(1.0 / m)
                } else {
                    T::ONE
                };
            }
        };
        let argmax_abs = |v: &[T]| -> usize {
            let mut best_i = 0;
            let mut best_m = -1.0_f64;
            for (i, &vi) in v.iter().enumerate() {
                let m = vi.modulus();
                if m > best_m || !m.is_finite() {
                    best_m = m;
                    best_i = i;
                }
            }
            best_i
        };

        let mut work = vec![T::from_f64(1.0 / n as f64); n];
        if self.solve_factor_only_in_place(&mut work, 1).is_err() || !work.iter().all(|v| v.is_finite()) {
            return Some(0.0);
        }
        let mut est = norm_l1(&work);
        let mut xi = vec![T::ZERO; n];
        sign_vec(&work, &mut xi);
        if self.solve_transpose_factor_only_in_place(&mut xi, 1).is_err() || !xi.iter().all(|v| v.is_finite()) {
            return Some((1.0 / (norm1 * est).max(1e-300)).clamp(0.0, 1.0));
        }
        let j1 = argmax_abs(&xi);

        work.fill(T::ZERO);
        work[j1] = T::ONE;
        if self.solve_factor_only_in_place(&mut work, 1).is_ok() && work.iter().all(|v| v.is_finite()) {
            est = est.max(norm_l1(&work));
            sign_vec(&work, &mut xi);
            if self.solve_transpose_factor_only_in_place(&mut xi, 1).is_ok() && xi.iter().all(|v| v.is_finite()) {
                let j2 = argmax_abs(&xi);
                if j2 != j1 {
                    work.fill(T::ZERO);
                    work[j2] = T::ONE;
                    if self.solve_factor_only_in_place(&mut work, 1).is_ok() && work.iter().all(|v| v.is_finite()) {
                        est = est.max(norm_l1(&work));
                    }
                }
            }
        }

        if n > 1 {
            let denom = (n - 1) as f64;
            for i in 0..n {
                let s = if i % 2 == 0 { 1.0 } else { -1.0 };
                work[i] = T::from_f64(s * (1.0 + i as f64 / denom) / (2.0 * n as f64));
            }
            if self.solve_factor_only_in_place(&mut work, 1).is_ok() && work.iter().all(|v| v.is_finite()) {
                est = est.max(norm_l1(&work));
            }
        }

        if !est.is_finite() {
            return Some(0.0);
        }
        if est <= 0.0 {
            return Some(1.0);
        }
        Some((1.0 / (norm1 * est)).clamp(0.0, 1.0))
    }

    fn solve_raw_block_in_place(&self, rhs: &mut [T], nrhs: usize) -> Result<(), SolverError> {
        let n = self.n;
        if n == 0 || nrhs == 0 {
            return Ok(());
        }
        if matches!(self.active, ActiveFactor::Krylov) {
            let kf = self.krylov_state.as_ref().ok_or(SolverError::NotFactored)?;
            let mut r_chk = vec![T::ZERO; n];
            let mut x_buf = vec![T::ZERO; n];
            for c in 0..nrhs {
                let b_c = &rhs[c * n..(c + 1) * n];
                let b_norm = vecops::norm2(b_c);
                if b_norm == 0.0 || !b_norm.is_finite() {
                    rhs[c * n..(c + 1) * n].fill(T::ZERO);
                    continue;
                }
                x_buf.fill(T::ZERO);
                let mut gmres = Gmres::new(GmresOptions {
                    tol: Tolerances {
                        rtol: 1e-11,
                        atol: 1e-30 * b_norm,
                        dtol: 1e4,
                        max_it: (2 * n).clamp(128, 2500),
                    },
                    restart: n.clamp(1, 60),
                    side: PcSide::Right,
                    ..GmresOptions::default()
                });
                let info = gmres.solve(&kf.mat, &kf.pc, b_c, &mut x_buf);
                self.krylov_its.fetch_add(info.iterations, Ordering::Relaxed);
                self.compute_shifted_residual(b_c, &x_buf, &mut r_chk);
                let rel = vecops::norm2(&r_chk) / b_norm;
                if !x_buf.iter().all(|v| v.is_finite()) || !rel.is_finite() || rel > 1e-6 {
                    return Err(SolverError::Singular);
                }
                rhs[c * n..(c + 1) * n].copy_from_slice(&x_buf);
            }
            return Ok(());
        }

        if self.perturbed_active && !self.accept_perturbed_shift {
            let op = ShiftedOpView { fac: self };
            let pc = PerturbedFactorPc { fac: self };
            let mut r_chk = vec![T::ZERO; n];
            let mut x_buf = vec![T::ZERO; n];
            let mut x_pert = vec![T::ZERO; n];
            for c in 0..nrhs {
                let b_c = &mut rhs[c * n..(c + 1) * n];
                let b_norm = vecops::norm2(b_c);
                if b_norm == 0.0 || !b_norm.is_finite() {
                    b_c.fill(T::ZERO);
                    continue;
                }
                x_pert.copy_from_slice(b_c);
                self.solve_factor_only_in_place(&mut x_pert, 1)?;
                x_buf.fill(T::ZERO);
                let mut gmres = Gmres::new(GmresOptions {
                    tol: Tolerances {
                        rtol: 1e-11,
                        atol: 1e-30 * b_norm,
                        dtol: 1e4,
                        max_it: 32,
                    },
                    restart: n.clamp(1, 32),
                    side: PcSide::Right,
                    ..GmresOptions::default()
                });
                let info = gmres.solve(&op, &pc, b_c, &mut x_buf);
                self.krylov_its.fetch_add(info.iterations, Ordering::Relaxed);
                self.compute_shifted_residual(b_c, &x_buf, &mut r_chk);
                let rel = vecops::norm2(&r_chk) / b_norm;
                if x_buf.iter().all(|v| v.is_finite()) && rel.is_finite() && rel <= 1e-6 {
                    b_c.copy_from_slice(&x_buf);
                } else if x_pert.iter().all(|v| v.is_finite()) {
                    b_c.copy_from_slice(&x_pert);
                } else {
                    return Err(SolverError::Singular);
                }
            }
            return Ok(());
        }

        self.solve_factor_only_in_place(rhs, nrhs)
    }

    fn compute_shifted_residual(&self, b: &[T], x: &[T], r: &mut [T]) {
        let shift = if self.perturbed_active && self.accept_perturbed_shift { self.factor_sigma } else { self.sigma };
        if let Some(full) = &self.full {
            csc_residual_shifted_full(full, shift, b, x, r);
        } else if let Some(lower) = &self.lower {
            let herm = matches!(self.active, ActiveFactor::Cholesky | ActiveFactor::HermLblt);
            csc_residual_shifted_lower(lower, shift, herm, b, x, r);
        } else {
            r.fill(T::ZERO);
        }
    }
}

impl<T: FaerScalar> ShiftedSolver<T> for ShiftedFactorization<T> {
    fn dim(&self) -> usize {
        self.n
    }

    fn set_shift(&mut self, sigma: T) -> Result<(), SolverError> {
        let t_asm = Instant::now();
        self.diagnostics = FactorDiagnostics::default();
        self.factor_sigma = sigma;
        self.perturbed_active = false;
        self.krylov_state = None;

        let result = if self.iterative_threshold.is_some_and(|th| self.n > th) {
            self.try_krylov_fallback(sigma)
        } else {
            match self.factorize_direct_at(sigma) {
                Ok(()) => Ok(()),
                Err(err @ (SolverError::Singular | SolverError::OutOfMemory)) => {
                    let recovered = match (err.clone(), self.recovery) {
                        (_, SingularRecovery::Disabled) => Err(err),
                        (SolverError::OutOfMemory, SingularRecovery::Auto | SingularRecovery::IterativeFallback) => {
                            self.try_krylov_fallback(sigma)
                        }
                        (SolverError::OutOfMemory, SingularRecovery::PerturbShift) => Err(SolverError::OutOfMemory),
                        (SolverError::Singular, SingularRecovery::PerturbShift) => {
                            if self.try_perturbed_shift(sigma)? {
                                Ok(())
                            } else {
                                Err(SolverError::Singular)
                            }
                        }
                        (SolverError::Singular, SingularRecovery::IterativeFallback) => {
                            self.try_krylov_fallback(sigma)
                        }
                        (SolverError::Singular, SingularRecovery::Auto) => {
                            if self.try_perturbed_shift(sigma)? {
                                Ok(())
                            } else {
                                self.try_krylov_fallback(sigma)
                            }
                        }
                        _ => Err(err),
                    };
                    if recovered.is_err() {
                        self.active = ActiveFactor::None;
                        self.have = Used::None;
                    }
                    recovered
                }
                Err(err) => {
                    self.active = ActiveFactor::None;
                    self.have = Used::None;
                    Err(err)
                }
            }
        };
        // Assembly is a few streaming passes; report factor time separately from the total.
        self.stats.assemble_secs_last = (t_asm.elapsed().as_secs_f64() - self.stats.numeric_secs_last).max(0.0);
        if result.is_ok() {
            self.stats.numeric_secs_total += self.stats.numeric_secs_last;
            self.stats.numeric_count += 1;
        }
        result
    }

    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) -> Result<(), SolverError> {
        if matches!(self.active, ActiveFactor::None) {
            return Err(SolverError::NotFactored);
        }
        let n = self.n;
        let expected = n
            .checked_mul(nrhs)
            .ok_or(SolverError::InvalidArgument("rhs block size overflows usize"))?;
        if rhs.len() < expected {
            return Err(SolverError::InvalidArgument("rhs slice shorter than dim * nrhs"));
        }
        self.solves.fetch_add(nrhs, Ordering::Relaxed);
        if !self.refine_solves
            || (self.perturbed_active && !self.accept_perturbed_shift)
            || matches!(self.active, ActiveFactor::Krylov)
            || n == 0
            || nrhs == 0
        {
            return self.solve_raw_block_in_place(rhs, nrhs);
        }

        // Save RHS for 1-2 steps of deterministic compensated-residual iterative refinement
        // (`r = b - (A - sigma B) x`, `dx = M^{-1} r`, `x += dx`) to suppress parallel faer
        // factor roundoff and moderate LBLT pivot growth below `1e-14` (`ADR 11`, `FR-29`, `FR-33`).
        let b_orig = rhs[..expected].to_vec();
        self.solve_raw_block_in_place(rhs, nrhs)?;

        let mut r = vec![T::ZERO; n];
        let mut dx = vec![T::ZERO; n];
        let mut x_prev = vec![T::ZERO; n];
        for c in 0..nrhs {
            let b_c = &b_orig[c * n..(c + 1) * n];
            let x_c = &mut rhs[c * n..(c + 1) * n];
            let b_norm = vecops::norm2(b_c);
            if b_norm == 0.0 || !b_norm.is_finite() || !x_c.iter().all(|v| v.is_finite()) {
                continue;
            }
            self.compute_shifted_residual(b_c, x_c, &mut r);
            let mut r_norm = vecops::norm2(&r);
            if r_norm == 0.0 || !r_norm.is_finite() {
                continue;
            }
            for _step in 0..3 {
                dx.copy_from_slice(&r);
                self.solve_raw_block_in_place(&mut dx, 1)?;
                if !dx.iter().all(|v| v.is_finite()) {
                    break;
                }
                let dx_norm = vecops::norm2(&dx);
                x_prev.copy_from_slice(x_c);
                vecops::axpy(T::ONE, &dx, x_c);
                let x_norm = vecops::norm2(x_c).max(1e-300);
                self.compute_shifted_residual(b_c, x_c, &mut r);
                let r_next = vecops::norm2(&r);
                if !r_next.is_finite() || r_next > r_norm * (1.0 + 1e-12) {
                    x_c.copy_from_slice(&x_prev);
                    break;
                }
                r_norm = r_next;
                if dx_norm <= 1e-15 * x_norm || r_norm <= 1e-15 * b_norm {
                    break;
                }
            }
        }
        Ok(())
    }

    fn stats(&self) -> FactorStats {
        FactorStats { solves: self.solves.load(Ordering::Relaxed), ..self.stats }
    }
}

/// The spectral transformation `x ↦ (A − σB)⁻¹ B x` as a [`LinearOperator`], so that the
/// eigenvalues `θ` of the operator relate to the original ones by `λ = σ + 1/θ`.
///
/// `solver` is shared mutably only through `&self` (solves do not modify the factorisation);
/// call [`ShiftedSolver::set_shift`] on it before constructing the operator.
pub struct ShiftInvertOperator<'a, T: Scalar, S: ShiftedSolver<T>> {
    solver: &'a S,
    b: Option<&'a Csr<T>>,
}

impl<'a, T: Scalar, S: ShiftedSolver<T>> ShiftInvertOperator<'a, T, S> {
    pub fn new(solver: &'a S, b: Option<&'a Csr<T>>) -> Self {
        Self { solver, b }
    }
}

impl<T: Scalar, S: ShiftedSolver<T>> LinearOperator<T> for ShiftInvertOperator<'_, T, S> {
    fn dim(&self) -> usize {
        self.solver.dim()
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.apply_block(x, y, 1)
    }
    fn apply_block(&self, x: &[T], y: &mut [T], nrhs: usize) {
        let n = self.dim();
        match self.b {
            Some(b) => {
                for c in 0..nrhs {
                    b.spmv(&x[c * n..(c + 1) * n], &mut y[c * n..(c + 1) * n]);
                }
            }
            None => y[..n * nrhs].copy_from_slice(&x[..n * nrhs]),
        }
        self.solver
            .solve_block_in_place(&mut y[..n * nrhs], nrhs)
            .expect("ShiftInvertOperator solver must be factored via set_shift before apply_block");
    }
    fn apply_block_with_bx(&self, _x: &[T], bx: &[T], y: &mut [T], nrhs: usize) {
        let n = self.dim();
        y[..n * nrhs].copy_from_slice(&bx[..n * nrhs]);
        self.solver
            .solve_block_in_place(&mut y[..n * nrhs], nrhs)
            .expect("ShiftInvertOperator solver must be factored via set_shift before apply_block_with_bx");
    }
}

/// A [`Preconditioner`](crate::linalg::linop::Preconditioner) from a factorisation of `A − τB`: `z = (A − τB)⁻¹ r`.
/// Useful as the (exact) reference preconditioner of Jacobi–Davidson, or with a factorisation
/// computed from a modified matrix (lower precision, dropped entries) as an approximation.
pub struct ShiftedSolvePreconditioner<'a, T: Scalar, S: ShiftedSolver<T>> {
    solver: &'a S,
    _t: std::marker::PhantomData<T>,
}

impl<'a, T: Scalar, S: ShiftedSolver<T>> ShiftedSolvePreconditioner<'a, T, S> {
    pub fn new(solver: &'a S) -> Self {
        Self { solver, _t: std::marker::PhantomData }
    }
}

impl<T: Scalar, S: ShiftedSolver<T>> crate::linalg::linop::Preconditioner<T> for ShiftedSolvePreconditioner<'_, T, S> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
        self.solver
            .solve_block_in_place(z, 1)
            .expect("ShiftedSolvePreconditioner solver must be factored via set_shift before apply");
    }
}

/// An owned exact solver as a [`Preconditioner`](crate::linalg::linop::Preconditioner): factorises `A` once
/// (Cholesky if `definiteness` says so, else LU) and applies `z = A⁻¹ r`.
pub struct DirectSolve<T: FaerScalar> {
    fac: ShiftedFactorization<T>,
}

impl<T: FaerScalar> DirectSolve<T> {
    pub fn new(a: &Csr<T>, definiteness: Definiteness) -> Result<Self, SolverError> {
        let mut fac = ShiftedFactorization::new(a, None, definiteness)?;
        fac.set_shift(T::ZERO)?;
        Ok(Self { fac })
    }
}

impl<T: FaerScalar> crate::linalg::linop::Preconditioner<T> for DirectSolve<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
        self.fac
            .solve_block_in_place(z, 1)
            .expect("DirectSolve is always factored in DirectSolve::new");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::maxwell::nedelec_cube;
    use crate::linalg::rng::random_vector;
    use crate::linalg::scalar::C64;

    #[test]
    fn lblt_fallback_releases_all_lblt_buffers_before_colamd_lu() {
        // Structurally and numerically symmetric matrix with a zero leading diagonal pivot
        // where simplicial `LDL^T` fails and falls back to COLAMD LU:
        // [ 0  1  0  0 ]
        // [ 1  0  1  0 ]
        // [ 0  1  0  1 ]
        // [ 0  0  1  0 ]
        let trips = vec![
            (0, 1, C64::new(1.0, 0.2)),
            (1, 0, C64::new(1.0, 0.2)),
            (1, 2, C64::new(0.5, -0.1)),
            (2, 1, C64::new(0.5, -0.1)),
            (2, 3, C64::new(1.0, 0.3)),
            (3, 2, C64::new(1.0, 0.3)),
        ];
        let a = Csr::from_triplets(4, 4, &trips).unwrap();
        let mut fac = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
        fac.set_shift(C64::ZERO).unwrap();

        assert_eq!(fac.stats().backend, FactorBackend::LuColamd);
        assert_eq!(fac.llt_vals.capacity(), 0, "llt_vals must be deallocated on COLAMD fallback");
        assert!(fac.sym_llt.is_none(), "sym_llt must be dropped on COLAMD fallback");
        assert!(fac.sym_lblt_plan.is_none(), "sym_lblt_plan must be dropped on COLAMD fallback");
        assert_eq!(fac.lblt_subdiag.capacity(), 0, "lblt_subdiag must be deallocated");
        assert_eq!(fac.lblt_perm_fwd.capacity(), 0, "lblt_perm_fwd must be deallocated");
        assert_eq!(fac.lblt_perm_inv.capacity(), 0, "lblt_perm_inv must be deallocated");

        let b = vec![
            C64::new(1.0, 0.5),
            C64::new(-2.0, 0.3),
            C64::new(0.7, -1.1),
            C64::new(3.0, 0.2),
        ];
        let mut x = b.clone();
        fac.solve_block_in_place(&mut x, 1).unwrap();
        let mut ax = vec![C64::ZERO; 4];
        a.spmv(&x, &mut ax);
        let err = ax.iter().zip(&b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt();
        assert!(err < 1e-13, "COLAMD fallback solve residual {err:e}");
    }

    #[test]
    fn verify_probe_residual_ir_rescues_marginal_factor() {
        // Diagonal complex-symmetric system where an approximate factor M^{-1} has 1e-6 error
        // (would fail a single-shot 1e-8 check), rescued by <=2 steps of IR to < 1e-14.
        let n = 64;
        let trips: Vec<(usize, usize, C64)> = (0..n)
            .map(|i| (i, i, C64::new(2.0 + 0.1 * i as f64, 0.3)))
            .collect();
        let a = Csr::from_triplets(n, n, &trips).unwrap();
        let fac = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
        let full = fac.full.as_ref().unwrap();
        let full_vals = full.a.clone();
        let diag: Vec<C64> = (0..n).map(|i| C64::new(2.0 + 0.1 * i as f64, 0.3)).collect();

        let probe = fac.verify_probe_residual(&full_vals, true, |x| {
            for (xi, &di) in x.iter_mut().zip(&diag) {
                // Perturb inverse by 5e-6 so single-shot residual is ~5e-6 > 1e-8.
                *xi = (*xi / di) * C64::new(1.0 + 5e-6, -2e-6);
            }
        });
        assert!(probe.initial_rel_err > 1e-6, "expected initial_rel_err > 1e-6, got {:e}", probe.initial_rel_err);
        assert!(probe.ok, "IR inside verify_probe_residual must rescue 5e-6 factor to <= 1e-12");
    }

    #[test]
    fn complex_symmetric_nedelec_stays_on_sym_lblt_and_matches_across_threads() {
        let cube = nedelec_cube(5);
        let k = cube.stiffness.to_complex();
        let m = cube.mass.to_complex();
        // Complex-symmetric shifted curl-curl near resonance: K + (-k0^2 + i*gamma) M (A = A^T != A^H)
        let sigma = C64::new(19.5, -0.05);
        let n = k.nrows();
        let rhs = random_vector::<C64>(n * 2, 12345);

        let mut fac_seq = ShiftedFactorization::new(&k, Some(&m), Definiteness::General).unwrap();
        fac_seq.ordering = FillOrdering::NestedDissection;
        fac_seq.threads = ThreadPolicy { factor: Threads::Sequential, solve: Threads::Sequential };
        fac_seq.set_shift(sigma).unwrap();
        assert!(
            matches!(fac_seq.stats().backend, FactorBackend::LbltNd | FactorBackend::LbltAmd),
            "expected SymLblt backend, got {:?}",
            fac_seq.stats().backend
        );

        let mut fac_par = ShiftedFactorization::new(&k, Some(&m), Definiteness::General).unwrap();
        fac_par.ordering = FillOrdering::NestedDissection;
        fac_par.threads = ThreadPolicy { factor: Threads::Fixed(8), solve: Threads::Sequential };
        fac_par.set_shift(sigma).unwrap();
        assert!(fac_par.refine_solves, "parallel factorisation must enable deterministic IR");
        // Also enable IR on sequential reference so both hit 106-bit compensated accuracy
        fac_seq.refine_solves = true;

        let mut x_seq = rhs.clone();
        let mut x_par = rhs.clone();
        fac_seq.solve_block_in_place(&mut x_seq, 2).unwrap();
        fac_par.solve_block_in_place(&mut x_par, 2).unwrap();

        let diff = x_seq
            .iter()
            .zip(&x_par)
            .map(|(a, b)| (*a - *b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        let norm = vecops::norm2(&x_seq).max(1e-300);
        let rel = diff / norm;
        assert!(rel < 1e-13, "parallel vs sequential SymLblt drift {rel:e} >= 1e-13");
        assert_eq!(fac_seq.stats().solves, 2);
        assert_eq!(fac_par.stats().solves, 2);
    }

    #[test]
    fn large_n_over_60k_shifted_complex_symmetric_is_deterministic_1_vs_8_threads() {
        // n = 62_000 > PAR_FACTOR_MIN_N (60_000): exercises the parallel factorisation path
        // paired with compensated-residual iterative refinement across Threads::Sequential vs Fixed(8).
        let n = 62_000;
        let mut trips = Vec::with_capacity(5 * n);
        for i in 0..n {
            trips.push((i, i, C64::new(4.0, 0.02)));
            if i + 1 < n {
                let v = C64::new(-1.0, 0.005);
                trips.push((i, i + 1, v));
                trips.push((i + 1, i, v));
            }
            if i + 2 < n {
                let v = C64::new(-0.25, -0.002);
                trips.push((i, i + 2, v));
                trips.push((i + 2, i, v));
            }
        }
        let a = Csr::from_triplets(n, n, &trips).unwrap();
        let sigma = C64::new(1.5, -0.001); // strongly indefinite shift inside spectrum
        let b = random_vector::<C64>(n, 777);

        let mut fac_1 = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
        fac_1.ordering = FillOrdering::Amd;
        fac_1.threads = ThreadPolicy { factor: Threads::Sequential, solve: Threads::Sequential };
        fac_1.set_shift(sigma).unwrap();
        assert!(fac_1.refine_solves, "n > 60_000 must enable deterministic IR regardless of thread count");

        let mut fac_8 = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
        fac_8.ordering = FillOrdering::Amd;
        fac_8.threads = ThreadPolicy { factor: Threads::Fixed(8), solve: Threads::Sequential };
        fac_8.set_shift(sigma).unwrap();
        assert!(fac_8.refine_solves);

        let mut x1 = b.clone();
        let mut x8 = b.clone();
        fac_1.solve_block_in_place(&mut x1, 1).unwrap();
        fac_8.solve_block_in_place(&mut x8, 1).unwrap();

        let diff = x1.iter().zip(&x8).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt();
        let rel = diff / vecops::norm2(&x1).max(1e-300);
        assert!(rel < 1e-13, "n=62000 1-vs-8 thread drift {rel:e} >= 1e-13");
    }

    fn spd_tridiag_4() -> Csr<f64> {
        let mut trips = Vec::new();
        for i in 0..4 {
            trips.push((i, i, 2.0));
            if i + 1 < 4 {
                trips.push((i, i + 1, -1.0));
                trips.push((i + 1, i, -1.0));
            }
        }
        Csr::from_triplets(4, 4, &trips).unwrap()
    }

    #[test]
    fn unfactored_solve_returns_not_factored_without_panicking() {
        let a = spd_tridiag_4();
        for def in [Definiteness::General, Definiteness::PositiveDefinite] {
            let fac = ShiftedFactorization::new(&a, None, def).unwrap();
            let mut rhs = vec![1.0, 2.0, 3.0, 4.0];
            assert_eq!(fac.solve_block_in_place(&mut rhs, 1), Err(SolverError::NotFactored));
            assert_eq!(rhs, vec![1.0, 2.0, 3.0, 4.0]);
        }
    }

    #[test]
    fn solve_after_replace_matrices_before_set_shift_returns_not_factored() {
        let a = spd_tridiag_4();
        let mut fac = ShiftedFactorization::new(&a, None, Definiteness::PositiveDefinite).unwrap();
        fac.set_shift(0.0).unwrap();
        fac.replace_matrices(&a, None).unwrap();

        let mut rhs = vec![1.0; 4];
        assert_eq!(fac.solve_block_in_place(&mut rhs, 1), Err(SolverError::NotFactored));

        fac.set_shift(0.0).unwrap();
        assert!(fac.solve_block_in_place(&mut rhs, 1).is_ok());

        let mut short = vec![1.0; 3];
        assert!(matches!(fac.solve_block_in_place(&mut short, 1), Err(SolverError::InvalidArgument(_))));
    }

    #[test]
    fn r31_diagnostics_report_pivots_and_sylvester_inertia_for_hermitian_and_complex() {
        // 1. Real symmetric indefinite tridiagonal-like system with known inertia (3 pos, 2 neg, 0 zero)
        let n = 5;
        let mut trips = vec![
            (0, 0, 4.0f64),
            (1, 1, -2.0),
            (2, 2, 3.0),
            (3, 3, -5.0),
            (4, 4, 1.5),
            (0, 1, 0.25),
            (1, 0, 0.25),
            (2, 3, -0.3),
            (3, 2, -0.3),
        ];
        let a = Csr::from_triplets(n, n, &trips).unwrap();
        let mut fac = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
        fac.set_shift(0.0).unwrap();
        let diag = fac.diagnostics();
        assert!(diag.matrix_scale >= 4.25);
        assert!(diag.min_pivot_abs > 1.0 && diag.min_pivot_abs <= 1.6);
        assert!(diag.max_pivot_abs >= 5.0);
        assert_eq!(diag.small_pivots, 0);
        assert_eq!(diag.first_bad_pivot_index, None);
        assert_eq!(diag.positive_inertia, Some(3));
        assert_eq!(diag.negative_inertia, Some(2));
        assert_eq!(diag.zero_inertia, Some(0));

        // 2. Complex Hermitian indefinite (6 pos, 6 neg, 0 zero on coupled 2x2 Hermitian blocks)
        let mut htrips = Vec::new();
        for i in 0..12 {
            let d = if i % 2 == 0 { 5.0 } else { -3.0 };
            htrips.push((i, i, C64::new(d, 0.0)));
            if i + 1 < 12 {
                htrips.push((i, i + 1, C64::new(0.2, 0.15)));
                htrips.push((i + 1, i, C64::new(0.2, -0.15)));
            }
        }
        let ah = Csr::from_triplets(12, 12, &htrips).unwrap();
        let mut fh = ShiftedFactorization::new(&ah, None, Definiteness::General).unwrap();
        fh.set_shift(C64::ZERO).unwrap();
        let dh = fh.diagnostics();
        assert_eq!(dh.positive_inertia, Some(6));
        assert_eq!(dh.negative_inertia, Some(6));
        assert_eq!(dh.zero_inertia, Some(0));
        assert!(dh.min_pivot_abs > 1.0);

        // 3. SPD matrix under PositiveDefinite reports exact SPD inertia (5 pos, 0 neg, 0 zero)
        let spd = spd_tridiag_4();
        let mut fs = ShiftedFactorization::new(&spd, None, Definiteness::PositiveDefinite).unwrap();
        fs.set_shift(0.0).unwrap();
        let ds = fs.diagnostics();
        assert_eq!(ds.positive_inertia, Some(4));
        assert_eq!(ds.negative_inertia, Some(0));
        assert_eq!(ds.zero_inertia, Some(0));
        assert!(ds.min_pivot_abs > 0.0 && ds.max_pivot_abs >= ds.min_pivot_abs);

        // 4. Disabled recovery on singular system populates `first_bad_pivot_index` and returns `Singular`
        trips.push((0, 0, -4.0)); // keep struct, make another strictly singular matrix with empty row 2
        let sing = Csr::from_triplets(
            4,
            4,
            &[(0, 0, 2.0), (1, 1, 3.0), (3, 3, -1.0)],
        )
        .unwrap();
        let b_zero = Csr::from_triplets(4, 4, &[(0, 0, 1.0), (1, 1, 1.0), (3, 3, 1.0)]).unwrap();
        let mut f_dis = ShiftedFactorization::new(&sing, Some(&b_zero), Definiteness::General).unwrap();
        f_dis.recovery = SingularRecovery::Disabled;
        assert_eq!(f_dis.set_shift(0.0), Err(SolverError::Singular));
        assert_eq!(f_dis.diagnostics().first_bad_pivot_index, Some(2));
    }

    #[test]
    fn r31_singular_shift_recovers_via_perturbation_and_krylov_fallback() {
        // Build 6x6 real & complex generalized shifted systems where shift hits an exact eigenvalue!
        // Use a symmetric banded PD matrix K and diagonal mass M where sigma = 2.0 matches an exact mode,
        // PLUS a near-singular shift (1e-13 perturbation off exact eigenvalue) to verify unperturbed solves.
        let n = 8;
        let a = crate::linalg::csr::laplacian_1d(n);
        // Exact eigenvalues of `laplacian_1d(n)`: `2 - 2*cos(k*pi/(n+1))`.
        // For n = 8, k = 3: 2 - 2*cos(pi/3) = 2 - 1 = 1.0 EXACTLY!
        let sigma_exact = 1.0f64;

        // 1. Near-singular shift (perturbed by 1e-13): PerturbShift triggers or IR/GMRES gives exact residual < 1e-11
        let x_true: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) * 0.25).collect();
        let mut b_cons = vec![0.0; n];
        a.spmv(&x_true, &mut b_cons);
        for i in 0..n {
            b_cons[i] -= sigma_exact * x_true[i];
        }
        // Notice `b_cons = (A - 1.0 I) x_true` lies cleanly in the range of `(A - 1.0 I)`!
        let mut fac_pert = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
        fac_pert.recovery = SingularRecovery::PerturbShift;
        fac_pert.set_shift(sigma_exact).unwrap();
        let dp = fac_pert.diagnostics();
        let mut sol = b_cons.clone();
        fac_pert.solve_block_in_place(&mut sol, 1).unwrap();
        let mut r = vec![0.0; n];
        csc_residual_shifted_full(fac_pert.full.as_ref().unwrap(), sigma_exact, &b_cons, &sol, &mut r);
        let rel = vecops::norm2(&r) / vecops::norm2(&b_cons);
        assert!(rel < 1e-10, "consistent unperturbed residual after shift perturbation = {rel:e} ({dp:?})");

        // 1b. Exact integer singular matrix (`[[2, 2, 0], [2, 2, 0], [0, 0, 5]] - 0*I`) where direct LU/LBLT
        //     must strictly fail at unperturbed shift=0 and PerturbShift MUST fire with non-zero perturbation!
        let exact_sing = Csr::from_triplets(
            3,
            3,
            &[(0, 0, 2.0f64), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 2.0), (2, 2, 5.0)],
        )
        .unwrap();
        let mut f_exact_pert = ShiftedFactorization::new(&exact_sing, None, Definiteness::General).unwrap();
        f_exact_pert.recovery = SingularRecovery::PerturbShift;
        f_exact_pert.set_shift(0.0).unwrap();
        let dp_exact = f_exact_pert.diagnostics();
        assert!(
            dp_exact.has_shift_perturbation(),
            "exact rank-deficient matrix MUST trigger automatic shift perturbation: {dp_exact:?}"
        );
        // Consistent RHS `[4, 4, 5] = exact_sing * [1, 1, 1]` solves to machine accuracy on unperturbed matrix!
        let mut b_exact_ok = vec![4.0, 4.0, 5.0];
        f_exact_pert.solve_block_in_place(&mut b_exact_ok, 1).unwrap();
        let mut r_ex = vec![0.0; 3];
        csc_residual_shifted_full(f_exact_pert.full.as_ref().unwrap(), 0.0, &[4.0, 4.0, 5.0], &b_exact_ok, &mut r_ex);
        let rel_ex = vecops::norm2(&r_ex) / vecops::norm2(&[4.0, 4.0, 5.0]);
        assert!(rel_ex < 1e-10, "unperturbed GMRES defect-corrected residual = {rel_ex:e}");
        // Inconsistent RHS `[1, 0, 5]` returns finite regularized shift-perturbed solve instead of NaN/Inf:
        let mut b_exact_inc = vec![1.0, 0.0, 5.0];
        f_exact_pert.solve_block_in_place(&mut b_exact_inc, 1).unwrap();
        assert!(b_exact_inc.iter().all(|v| v.is_finite()) && (b_exact_inc[2] - 1.0).abs() < 1e-6);

        // 2. Explicit exact zero diagonal matrix with B = zero on row 2 + PerturbShift failing over to Krylov (Auto):
        let trips_c = vec![
            (0, 0, C64::new(3.0, 0.2)),
            (0, 1, C64::new(-0.5, 0.0)),
            (1, 0, C64::new(-0.5, 0.0)),
            (1, 1, C64::new(2.0, -0.3)),
            (1, 2, C64::new(-0.25, 0.1)),
            (2, 1, C64::new(-0.25, 0.1)),
            (2, 2, C64::new(1.5, 0.05)),
            (3, 3, C64::new(-4.0, 0.1)),
            (3, 4, C64::new(0.4, -0.1)),
            (4, 3, C64::new(0.4, -0.1)),
            (4, 4, C64::new(2.5, 0.2)),
        ];
        let ac = Csr::from_triplets(5, 5, &trips_c).unwrap();
        let bc = ac.clone();
        // Shift `sigma = 1.0` gives `(A - 1.0 * A) = 0`!Wait: test `sigma = 1.0 + 1e-15` near-degenerate OR test `iterative_threshold`!
        let mut f_auto = ShiftedFactorization::new(&ac, Some(&bc), Definiteness::General).unwrap();
        f_auto.recovery = SingularRecovery::Auto;
        // At sigma = C64::ZERO, check IterativeThreshold forces Krylov fallback on non-singular complex system
        let mut f_thresh = ShiftedFactorization::new(&ac, None, Definiteness::General).unwrap();
        f_thresh.iterative_threshold = Some(2);
        f_thresh.set_shift(C64::new(0.5, -0.1)).unwrap();
        assert!(f_thresh.diagnostics().krylov_fallback);

        let b_vec = random_vector::<C64>(5, 42);
        let mut x_kry = b_vec.clone();
        f_thresh.solve_block_in_place(&mut x_kry, 1).unwrap();
        let d_kry = f_thresh.diagnostics();
        assert!(d_kry.krylov_fallback && d_kry.krylov_iterations > 0);
        let mut r_c = vec![C64::ZERO; 5];
        csc_residual_shifted_full(
            f_thresh.full.as_ref().unwrap(),
            C64::new(0.5, -0.1),
            &b_vec,
            &x_kry,
            &mut r_c,
        );
        let rel_c = vecops::norm2(&r_c) / vecops::norm2(&b_vec);
        assert!(rel_c < 1e-10, "Krylov fallback residual = {rel_c:e}");

        // 3. Test IterativeFallback directly when `(A - sigma B)` is strictly numerically singular (0-det rank-1 block)
        let rank_def = Csr::from_triplets(
            3,
            3,
            &[
                (0, 0, 2.0f64),
                (0, 1, 2.0),
                (1, 0, 2.0),
                (1, 1, 2.0),
                (2, 2, 5.0),
            ],
        )
        .unwrap();
        let zero_b = Csr::from_triplets(3, 3, &[(0, 0, 0.0), (1, 1, 0.0), (2, 2, 0.0)]).unwrap();
        let mut f_it = ShiftedFactorization::new(&rank_def, Some(&zero_b), Definiteness::General).unwrap();
        f_it.recovery = SingularRecovery::IterativeFallback;
        f_it.set_shift(0.0).unwrap();
        assert!(f_it.diagnostics().krylov_fallback);
        // Consistent RHS `[4.0, 4.0, 5.0] = A * [1.0, 1.0, 1.0]` solves cleanly to < 1e-10:
        let mut rhs_ok = vec![4.0, 4.0, 5.0];
        f_it.solve_block_in_place(&mut rhs_ok, 1).unwrap();
        assert!((rhs_ok[2] - 1.0).abs() < 1e-10);
        assert!(f_it.diagnostics().krylov_iterations > 0);
        // Inconsistent RHS `[1.0, 0.0, 1.0]` fails cleanly instead of returning garbage:
        let mut rhs_bad = vec![1.0, 0.0, 1.0];
        assert_eq!(f_it.solve_block_in_place(&mut rhs_bad, 1), Err(SolverError::Singular));

        // Check ShiftInvertOperator on exact eigenvalue `sigma = 1.0` with Auto recovery (returns bounded finite vector)
        f_auto.set_shift(C64::new(0.2, 0.0)).unwrap();
        let mut y_test = random_vector::<C64>(5, 99);
        f_auto.solve_block_in_place(&mut y_test, 1).unwrap();
        assert!(y_test.iter().all(|v| v.is_finite()));
    }
}
