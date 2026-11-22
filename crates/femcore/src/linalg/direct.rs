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
use crate::linalg::linop::LinearOperator;
use crate::linalg::ordering::{FillOrdering, adjacency_from_lower, nested_dissection};
use crate::linalg::scalar::Scalar;
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
use std::time::Instant;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolverError {
    /// The matrix is singular (or numerically so) at this shift.
    Singular,
    /// A Cholesky factorisation was required but the matrix is not positive definite.
    NotPositiveDefinite,
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

impl Threads {
    pub(crate) fn par(self) -> Par {
        match self {
            Threads::Global => {
                if crate::par::factor_threads() <= 1 {
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
    /// systems up to 60,000 unknowns factorise with `Par::Seq` so direct solves below the 3D AMG
    /// crossover threshold are 100% bit-deterministic across `NYX_THREADS` while avoiding supernode
    /// thread-synchronisation overhead on small/medium problems.
    pub(crate) fn par_for_factor(self, n: usize) -> Par {
        match self {
            Threads::Global if n <= 60_000 || crate::par::factor_threads() <= 1 => Par::Seq,
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
    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize);
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
        if last == Some((c, r)) {
            *av.last_mut().unwrap() += x;
            *bv.last_mut().unwrap() += y;
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

fn csc_spmv<T: Scalar>(n: usize, col_ptr: &[u32], row_idx: &[u32], vals: &[T], x: &[T], y: &mut [T]) {
    y.fill(T::ZERO);
    for c in 0..n {
        let xc = x[c];
        if xc == T::ZERO {
            continue;
        }
        let start = col_ptr[c] as usize;
        let end = col_ptr[c + 1] as usize;
        for k in start..end {
            y[row_idx[k] as usize] += vals[k] * xc;
        }
    }
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
) -> bool {
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
            return false;
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
                return false;
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
                return false;
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
    true
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
) -> bool {
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
        if !dense_sym_lblt(
            &mut dense_top[..s_ncols * s_ncols],
            s_ncols,
            &mut subdiag[s_start..s_end],
            &mut perm_fwd[s_start..s_end],
            &mut perm_inv[s_start..s_end],
            s_start,
        ) {
            return false;
        }
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
    true
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
) -> bool {
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
            return false;
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
    true
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
    vals: Vec<T>,
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
    have: Used,
    active: ActiveFactor,
    /// Thread policy for factorisation and solves (see [`ThreadPolicy`]).
    pub threads: ThreadPolicy,
    /// Fill-reducing ordering for Cholesky and symmetric indefinite `LBL^H`/`LBL^T`; set before the first factorisation.
    pub ordering: FillOrdering,
    stats: FactorStats,
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
        let nnz = full.as_ref().map_or_else(|| lower.as_ref().map_or(0, |l| 2 * l.row_idx.len()), |f| f.row_idx.len());
        Ok(Self {
            n,
            definiteness,
            full,
            tmap,
            lower,
            src,
            vals: vec![T::ZERO; nnz],
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
            have: Used::None,
            active: ActiveFactor::None,
            threads: ThreadPolicy::default(),
            ordering: FillOrdering::default(),
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
        })
    }

    /// Replaces the *values* of `A` (and `B`) without redoing the symbolic analysis, for matrices
    /// that change every step (Jacobians, `T(z) = Σ f_j(z) A_j`). The sparsity pattern of the union
    /// must be identical to the one the object was built with (otherwise `InvalidArgument`);
    /// [`ShiftedSolver::set_shift`] must be called again before solving.
    pub fn replace_matrices(&mut self, a: &Csr<T>, b: Option<&Csr<T>>) -> Result<(), SolverError> {
        self.ensure_full()?;
        let full = assemble(a, b, false)?;
        let old = self.full.as_ref().unwrap();
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
        Ok(())
    }

    /// Builds the full union pattern if only the lower triangle exists so far (needed by the LU path).
    fn ensure_full(&mut self) -> Result<(), SolverError> {
        if let Some(full) = &self.full {
            if self.tmap.is_none() {
                self.tmap = Some(build_transpose_map(full, self.n));
            }
        } else {
            let (a, b) = self.src.as_ref().expect("source matrices kept while the full pattern is missing");
            let f = assemble(a, b.as_ref(), false)?;
            let nnz = f.row_idx.len();
            self.tmap = Some(build_transpose_map(&f, self.n));
            self.full = Some(f);
            self.vals = vec![T::ZERO; nnz];
        }
        Ok(())
    }

    fn ensure_lower(&mut self) {
        if self.lower.is_none() {
            let full = self.full.as_ref().expect("full pattern must exist before extracting lower");
            self.lower = Some(extract_lower(full, self.n));
        }
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
        let lower = self.lower.as_ref().expect("lower pattern must exist");
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

    fn verify_probe_residual(&self, full_vals: &[T], solve_one: impl FnOnce(&mut [T])) -> bool {
        let n = self.n;
        let mut probe = Vec::with_capacity(n);
        let mut p_norm_sq = 0.0f64;
        for i in 0..n {
            let v = (((i.wrapping_mul(2654435761) ^ 0x5a5a_a5a5) % 10007) as f64 / 5003.5) - 1.0;
            let s = T::from_f64(v);
            p_norm_sq += s.modulus_sq();
            probe.push(s);
        }
        let p_norm = p_norm_sq.sqrt().max(1e-300);
        let mut x = probe.clone();
        solve_one(&mut x);
        if !x.iter().all(|v| v.is_finite()) {
            return false;
        }
        let full = self.full.as_ref().unwrap();
        let mut ax = vec![T::ZERO; n];
        csc_spmv(n, &full.col_ptr, &full.row_idx, full_vals, &x, &mut ax);
        let err_sq: f64 = ax.iter().zip(&probe).map(|(a, b)| (*a - *b).modulus_sq()).sum();
        let rel_err = err_sq.sqrt() / p_norm;
        rel_err.is_finite() && rel_err <= 1e-8
    }

    fn try_cholesky(&mut self, sigma: T) -> Result<bool, SolverError> {
        let n = self.n;
        let par = self.threads.factor.par_for_factor(n);
        self.ensure_sym_llt()?;
        let lower = self.lower.as_ref().expect("lower pattern exists for PositiveDefinite");
        let mut vals = std::mem::take(&mut self.vals);
        Self::shifted_values(lower, sigma, &mut vals);
        let sym = self.sym_llt.as_ref().unwrap();
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
        self.vals = vals;
        match res {
            Ok(_) => {
                self.have = Used::Cholesky;
                self.active = ActiveFactor::Cholesky;
                self.stats.used = Used::Cholesky;
                self.stats.backend = match self.sym_llt_ord {
                    FillOrdering::NestedDissection => FactorBackend::CholeskyNd,
                    _ => FactorBackend::CholeskyAmd,
                };
                self.stats.factor_nnz = self.sym_llt.as_ref().unwrap().len_val();
                Ok(true)
            }
            Err(_) => Ok(false),
        }
    }

    /// Supernodal/simplicial Hermitian `LBL^H` (`faer::intranode_lblt`) for real-symmetric or
    /// complex-Hermitian indefinite matrices, guarded by finiteness + 1-vector backward residual `<= 1e-8`.
    fn try_herm_lblt(&mut self, sigma: T, full_vals: &[T]) -> Result<bool, SolverError> {
        let n = self.n;
        if n == 0 {
            return Ok(false);
        }
        self.ensure_lower();
        self.ensure_sym_llt()?;
        let par = self.threads.factor.par_for_factor(n);
        let lower = self.lower.as_ref().unwrap();
        let mut lvals = Vec::with_capacity(lower.row_idx.len());
        Self::shifted_values(lower, sigma, &mut lvals);

        self.lblt_subdiag.resize(n, T::ZERO);
        self.lblt_perm_fwd.resize(n, 0);
        self.lblt_perm_inv.resize(n, 0);

        let sym = self.sym_llt.as_ref().unwrap();
        let mat = SparseColMatRef::new(Self::pattern(lower, n), T::faer_slice(&lvals));
        let t0 = Instant::now();
        let mut mem = MemBuffer::new(sym.factorize_numeric_intranode_lblt_scratch::<T::Elem>(par, Default::default()));
        let _lblt = sym.factorize_numeric_intranode_lblt::<T::Elem>(
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

        if !self.llt_vals.iter().all(|v| v.is_finite()) || !self.lblt_subdiag.iter().all(|v| v.is_finite()) {
            return Ok(false);
        }

        let solve_par = self.threads.solve.par();
        let ok = self.verify_probe_residual(full_vals, |x| {
            let sym = self.sym_llt.as_ref().unwrap();
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
        if !ok {
            return Ok(false);
        }

        self.stats.numeric_secs_last = fac_secs;
        self.have = Used::Lu;
        self.active = ActiveFactor::HermLblt;
        self.stats.used = Used::Lu;
        self.stats.backend = match self.sym_llt_ord {
            FillOrdering::NestedDissection => FactorBackend::LbltNd,
            _ => FactorBackend::LbltAmd,
        };
        self.stats.factor_nnz = self.sym_llt.as_ref().unwrap().len_val();
        Ok(true)
    }

    /// Non-conjugate supernodal/simplicial complex-symmetric `LBL^T` (`A = A^T != A^H`, e.g. PML /
    /// lossy dielectrics / surface impedance), reusing `faer`'s `SymbolicCholesky` (AMD/ND/Auto) and
    /// BLAS-3 dense kernels, guarded by finiteness + 1-vector backward residual `<= 1e-8`.
    fn try_sym_lblt(&mut self, sigma: T, full_vals: &[T]) -> Result<bool, SolverError> {
        let n = self.n;
        if n == 0 {
            return Ok(false);
        }
        self.ensure_lower();
        self.ensure_sym_llt()?;
        if self.sym_lblt_plan.is_none() {
            let t0 = Instant::now();
            let plan = build_sym_lblt_plan(self.lower.as_ref().unwrap(), self.sym_llt.as_ref().unwrap(), n);
            self.sym_lblt_plan = Some(plan);
            self.stats.analyse_secs += t0.elapsed().as_secs_f64();
        }
        let par = self.threads.factor.par_for_factor(n);
        let lower = self.lower.as_ref().unwrap();
        let mut lvals = Vec::with_capacity(lower.row_idx.len());
        Self::shifted_values(lower, sigma, &mut lvals);

        self.lblt_subdiag.resize(n, T::ZERO);
        self.lblt_perm_fwd.resize(n, 0);
        self.lblt_perm_inv.resize(n, 0);

        let t0 = Instant::now();
        let fac_ok = {
            let sym = self.sym_llt.as_ref().unwrap();
            let plan = self.sym_lblt_plan.as_ref().unwrap();
            match sym.raw() {
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
            }
        };
        let fac_secs = t0.elapsed().as_secs_f64();
        if !fac_ok || !self.llt_vals.iter().all(|v| v.is_finite()) || !self.lblt_subdiag.iter().all(|v| v.is_finite()) {
            return Ok(false);
        }

        let solve_par = self.threads.solve.par();
        let ok = self.verify_probe_residual(full_vals, |x| {
            let sym = self.sym_llt.as_ref().unwrap();
            match sym.raw() {
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
            }
        });
        if !ok {
            return Ok(false);
        }

        self.stats.numeric_secs_last = fac_secs;
        self.have = Used::Lu;
        self.active = ActiveFactor::SymLblt;
        self.stats.used = Used::Lu;
        self.stats.backend = match self.sym_llt_ord {
            FillOrdering::NestedDissection => FactorBackend::LbltNd,
            _ => FactorBackend::LbltAmd,
        };
        self.stats.factor_nnz = self.sym_llt.as_ref().unwrap().len_val();
        Ok(true)
    }

    fn do_colamd_lu(&mut self, full_vals: &[T]) -> Result<(), SolverError> {
        let n = self.n;
        let par = self.threads.factor.par_for_factor(n);
        if self.sym_lu.is_none() {
            let t0 = Instant::now();
            let full = self.full.as_ref().unwrap();
            let sym = factorize_symbolic_lu(Self::pattern(full, n), Default::default())
                .map_err(|_| SolverError::OutOfMemory)?;
            // Estimate COLAMD symmetric fill when structurally symmetric so factor_nnz is comparable.
            let est_nnz = if self.tmap.as_ref().is_some_and(|t| t.structurally_symmetric) {
                self.ensure_lower();
                let lower = self.lower.as_ref().unwrap();
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
            self.sym_lu = Some(sym);
            self.sym_lu_est_nnz = est_nnz;
            self.stats.analyse_secs += t0.elapsed().as_secs_f64();
            self.stats.symbolic_count += 1;
        }
        let sym = self.sym_lu.as_ref().unwrap();
        let mat = SparseColMatRef::new(Self::pattern(self.full.as_ref().unwrap(), n), T::faer_slice(full_vals));
        let t0 = Instant::now();
        let mut mem = MemBuffer::new(sym.factorize_numeric_lu_scratch::<T::Elem>(par, Default::default()));
        let res = sym.factorize_numeric_lu::<T::Elem>(
            &mut self.lu_num,
            mat,
            par,
            MemStack::new(&mut mem),
            Default::default(),
        );
        let ok = res.is_ok();
        clear_upper_vector_state();
        self.stats.numeric_secs_last = t0.elapsed().as_secs_f64();
        if ok {
            self.have = Used::Lu;
            self.active = ActiveFactor::ColamdLu;
            self.stats.used = Used::Lu;
            self.stats.backend = FactorBackend::LuColamd;
            self.stats.factor_nnz = self.sym_lu_est_nnz;
            Ok(())
        } else {
            Err(SolverError::Singular)
        }
    }

    fn do_lu(&mut self, sigma: T) -> Result<(), SolverError> {
        self.ensure_full()?;
        let mut vals = std::mem::take(&mut self.vals);
        Self::shifted_values(self.full.as_ref().unwrap(), sigma, &mut vals);

        let policy = backend_policy();
        let tmap = self.tmap.as_ref().unwrap();
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

        let res = match policy {
            BackendPolicy::ForceColamd => self.do_colamd_lu(&vals),
            BackendPolicy::ForceLblt => {
                if (is_herm && self.try_herm_lblt(sigma, &vals)?)
                    || (is_cplx_sym && self.try_sym_lblt(sigma, &vals)?)
                {
                    Ok(())
                } else {
                    self.do_colamd_lu(&vals)
                }
            }
            BackendPolicy::Auto => {
                // If COLAMD LU was already analysed for an earlier shift (e.g. unsymmetric or
                // previous fallback), keep reusing it so `symbolic_count` stays 1.
                if self.sym_lu.is_some() {
                    self.do_colamd_lu(&vals)
                } else if (is_herm && self.try_herm_lblt(sigma, &vals)?)
                    || (is_cplx_sym && self.try_sym_lblt(sigma, &vals)?)
                {
                    // Real symmetric / complex Hermitian (`LBL^H`) or complex-symmetric (`LBL^T`)
                    // supernodal Bunch-Kaufman with AMD/ND/Auto ordering.
                    Ok(())
                } else {
                    // Unsymmetric or LBL^T pivot fallback: faer COLAMD LU with full partial pivoting.
                    self.do_colamd_lu(&vals)
                }
            }
        };

        self.vals = vals;
        res
    }
}

impl<T: FaerScalar> ShiftedSolver<T> for ShiftedFactorization<T> {
    fn dim(&self) -> usize {
        self.n
    }

    fn set_shift(&mut self, sigma: T) -> Result<(), SolverError> {
        let t_asm = Instant::now();
        let mut done = false;
        let result = (|| {
            if self.definiteness == Definiteness::PositiveDefinite {
                if self.try_cholesky(sigma)? {
                    done = true;
                } else {
                    self.stats.cholesky_fallbacks += 1;
                }
            }
            if !done { self.do_lu(sigma) } else { Ok(()) }
        })();
        // Assembly is a few streaming passes; report factor time separately from the total.
        self.stats.assemble_secs_last = (t_asm.elapsed().as_secs_f64() - self.stats.numeric_secs_last).max(0.0);
        if result.is_ok() {
            self.stats.numeric_secs_total += self.stats.numeric_secs_last;
            self.stats.numeric_count += 1;
        }
        result
    }

    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) {
        let n = self.n;
        assert!(rhs.len() >= n * nrhs);
        let par = self.threads.solve.par();
        match self.active {
            ActiveFactor::Cholesky => {
                let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
                let sym = self.sym_llt.as_ref().unwrap();
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
                let sym = self.sym_llt.as_ref().unwrap();
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
                let sym = self.sym_llt.as_ref().unwrap();
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
                let sym = self.sym_lu.as_ref().unwrap();
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                LuRef::new_unchecked(sym, &self.lu_num).solve_in_place_with_conj(
                    Conj::No,
                    m,
                    par,
                    MemStack::new(&mut mem),
                );
            }
            ActiveFactor::None => panic!("set_shift must succeed before solving"),
        }
        clear_upper_vector_state();
    }

    fn stats(&self) -> FactorStats {
        self.stats
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
        self.solver.solve_block_in_place(&mut y[..n * nrhs], nrhs);
    }
    fn apply_block_with_bx(&self, _x: &[T], bx: &[T], y: &mut [T], nrhs: usize) {
        let n = self.dim();
        y[..n * nrhs].copy_from_slice(&bx[..n * nrhs]);
        self.solver.solve_block_in_place(&mut y[..n * nrhs], nrhs);
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
        self.solver.solve_block_in_place(z, 1);
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
        self.fac.solve_block_in_place(z, 1);
    }
}
