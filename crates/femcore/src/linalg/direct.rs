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
use faer::perm::PermRef;
use faer::sparse::linalg::cholesky::{LltRef, SymbolicCholesky, SymmetricOrdering, factorize_symbolic_cholesky};
use faer::sparse::linalg::lu::{LuRef, NumericLu, SymbolicLu, factorize_symbolic_lu};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use faer::{Conj, MatMut, Par, Side};
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
            Threads::Global => faer::get_global_parallelism(),
            Threads::Sequential => Par::Seq,
            Threads::Fixed(n) if n <= 1 => Par::Seq,
            Threads::Fixed(n) => Par::rayon(n),
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
    /// Cholesky attempts that failed and fell back to LU.
    pub cholesky_fallbacks: usize,
    /// Right-hand sides solved.
    pub solves: usize,
    /// Stored entries of the factor (`L` for Cholesky) once a Cholesky factorisation has been analysed; 0 before.
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

/// `faer`-backed [`ShiftedSolver`] with symbolic reuse.
pub struct ShiftedFactorization<T: FaerScalar> {
    n: usize,
    definiteness: Definiteness,
    /// Union pattern of the whole matrix; built on first use when the Cholesky path (which reads the lower triangle only) is enough.
    full: Option<Assembled<T>>,
    lower: Option<Assembled<T>>,
    /// Copies of the inputs kept only while `full` is not built (positive definite case).
    src: Option<(Csr<T>, Option<Csr<T>>)>,
    vals: Vec<T>,
    sym_lu: Option<SymbolicLu<u32>>,
    sym_llt: Option<SymbolicCholesky<u32>>,
    lu_num: NumericLu<u32, T::Elem>,
    llt_vals: Vec<T::Elem>,
    have: Used,
    /// Thread policy for factorisation and solves (see [`ThreadPolicy`]).
    pub threads: ThreadPolicy,
    /// Fill-reducing ordering of the Cholesky path; set before the first factorisation.
    pub ordering: FillOrdering,
    stats: FactorStats,
}

impl<T: FaerScalar> ShiftedFactorization<T> {
    /// Prepares the pattern for `A − σB` (`B = None` means the identity). With
    /// [`Definiteness::PositiveDefinite`], `A` and `B` must be symmetric/Hermitian (only their
    /// lower triangles are used by the Cholesky path).
    pub fn new(a: &Csr<T>, b: Option<&Csr<T>>, definiteness: Definiteness) -> Result<Self, SolverError> {
        let t0 = Instant::now();
        let pd = definiteness == Definiteness::PositiveDefinite;
        let (full, lower, src) = if pd {
            (None, Some(assemble_lower_sorted(a, b)?), Some((a.clone(), b.cloned())))
        } else {
            (Some(assemble(a, b, false)?), None, None)
        };
        let nnz = full.as_ref().map_or_else(|| lower.as_ref().map_or(0, |l| 2 * l.row_idx.len()), |f| f.row_idx.len());
        Ok(Self {
            n: a.nrows(),
            definiteness,
            full,
            lower,
            src,
            vals: vec![T::ZERO; nnz],
            sym_lu: None,
            sym_llt: None,
            lu_num: NumericLu::new(),
            llt_vals: Vec::new(),
            have: Used::None,
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
        if self.definiteness == Definiteness::PositiveDefinite {
            self.lower = Some(assemble_lower_sorted(a, b)?);
            self.src = Some((a.clone(), b.cloned()));
        }
        self.full = Some(full);
        self.have = Used::None;
        Ok(())
    }

    /// Builds the full union pattern if only the lower triangle exists so far (needed by the LU path).
    fn ensure_full(&mut self) -> Result<(), SolverError> {
        if self.full.is_none() {
            let (a, b) = self.src.as_ref().expect("source matrices kept while the full pattern is missing");
            self.full = Some(assemble(a, b.as_ref(), false)?);
            let nnz = self.full.as_ref().unwrap().row_idx.len();
            self.vals = vec![T::ZERO; nnz];
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

    fn try_cholesky(&mut self, sigma: T) -> Result<bool, SolverError> {
        let n = self.n;
        let par = self.threads.factor.par();
        let lower = self.lower.as_ref().expect("lower pattern exists for PositiveDefinite");
        if self.sym_llt.is_none() {
            let t0 = Instant::now();
            let pattern = Self::pattern(lower, n);
            let amd = || factorize_symbolic_cholesky(pattern, Side::Lower, SymmetricOrdering::Amd, Default::default());
            let nd = |perm: &(Vec<u32>, Vec<u32>)| {
                factorize_symbolic_cholesky(
                    pattern,
                    Side::Lower,
                    SymmetricOrdering::Custom(PermRef::new_checked(&perm.0, &perm.1, n)),
                    Default::default(),
                )
            };
            let make_perm = || {
                let adj = adjacency_from_lower(n, &lower.col_ptr, &lower.row_idx);
                let order = nested_dissection(&adj, 100);
                let fwd: Vec<u32> = order.iter().map(|&o| o as u32).collect();
                let mut inv = vec![0u32; n];
                for (new, &old) in order.iter().enumerate() {
                    inv[old] = new as u32;
                }
                (fwd, inv)
            };
            let sym = match self.ordering {
                FillOrdering::Amd => amd(),
                FillOrdering::NestedDissection => nd(&make_perm()),
                FillOrdering::Auto => {
                    // only worth trying for graphs with enough fill to matter
                    let a = amd().map_err(|_| SolverError::OutOfMemory)?;
                    if a.len_val() < 20 * n.max(1) {
                        Ok(a)
                    } else {
                        let b = nd(&make_perm()).map_err(|_| SolverError::OutOfMemory)?;
                        Ok(if b.len_val() < a.len_val() { b } else { a })
                    }
                }
            }
            .map_err(|_| SolverError::OutOfMemory)?;
            self.stats.factor_nnz = sym.len_val();
            self.llt_vals = vec![<T::Elem as faer::traits::ComplexField>::zero_impl(); sym.len_val()];
            self.sym_llt = Some(sym);
            self.stats.analyse_secs += t0.elapsed().as_secs_f64();
            self.stats.symbolic_count += 1;
        }
        let mut vals = std::mem::take(&mut self.vals);
        Self::shifted_values(lower, sigma, &mut vals);
        let sym = self.sym_llt.as_ref().unwrap();
        let mat = SparseColMatRef::new(Self::pattern(lower, n), T::faer_slice(&vals));
        let t0 = Instant::now();
        let mut mem = MemBuffer::new(sym.factorize_numeric_llt_scratch::<T::Elem>(par, Default::default()));
        let res = sym.factorize_numeric_llt::<T::Elem>(
            &mut self.llt_vals,
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
                self.stats.used = Used::Cholesky;
                Ok(true)
            }
            Err(_) => Ok(false),
        }
    }

    fn do_lu(&mut self, sigma: T) -> Result<(), SolverError> {
        self.ensure_full()?;
        let n = self.n;
        let par = self.threads.factor.par();
        if self.sym_lu.is_none() {
            let t0 = Instant::now();
            self.sym_lu = Some(
                factorize_symbolic_lu(Self::pattern(self.full.as_ref().unwrap(), n), Default::default())
                    .map_err(|_| SolverError::OutOfMemory)?,
            );
            self.stats.analyse_secs += t0.elapsed().as_secs_f64();
            self.stats.symbolic_count += 1;
        }
        let mut vals = std::mem::take(&mut self.vals);
        Self::shifted_values(self.full.as_ref().unwrap(), sigma, &mut vals);
        let sym = self.sym_lu.as_ref().unwrap();
        let mat = SparseColMatRef::new(Self::pattern(self.full.as_ref().unwrap(), n), T::faer_slice(&vals));
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
        self.vals = vals;
        if ok {
            self.have = Used::Lu;
            self.stats.used = Used::Lu;
            Ok(())
        } else {
            Err(SolverError::Singular)
        }
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
        self.stats.assemble_secs_last = t_asm.elapsed().as_secs_f64() - self.stats.numeric_secs_last;
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
        let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut rhs[..n * nrhs]), n, nrhs);
        match self.have {
            Used::Cholesky => {
                let sym = self.sym_llt.as_ref().unwrap();
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                LltRef::new(sym, &self.llt_vals).solve_in_place_with_conj(Conj::No, m, par, MemStack::new(&mut mem));
            }
            Used::Lu => {
                let sym = self.sym_lu.as_ref().unwrap();
                let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Elem>(nrhs, par));
                LuRef::new_unchecked(sym, &self.lu_num).solve_in_place_with_conj(
                    Conj::No,
                    m,
                    par,
                    MemStack::new(&mut mem),
                );
            }
            Used::None => panic!("set_shift must succeed before solving"),
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
