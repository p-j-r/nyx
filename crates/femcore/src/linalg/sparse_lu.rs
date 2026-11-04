//! A sparse LU factorisation that **accepts a user fill-reducing ordering** (port item S2, LU half).
//!
//! `faer`'s sparse LU chooses its own column ordering (COLAMD) inside the symbolic analysis and has no
//! hook for another one, whereas its Cholesky does accept a custom permutation (which is how the
//! nested-dissection ordering of [`crate::linalg::ordering`] is used, see `LEARNINGS.md` 24). That left the
//! indefinite and complex-*symmetric* matrices (PML, absorbing layers, shifted Helmholtz) without access to
//! the ordering that helps most on 3-D meshes. This module closes the gap with a small left-looking LU in the
//! style of Gilbert & Peierls (1988), written from the paper:
//!
//! * `P A Q = L U` with a column permutation `Q` chosen by the caller ([`LuOrdering`]): natural, faer's AMD
//!   on the pattern of `A + Aᵀ`, this crate's nested dissection, or any permutation.
//! * Column `k` of `L` and `U` comes from a *sparse* triangular solve whose non-zero pattern is found by a
//!   depth-first search on the graph of the columns of `L` already computed, so the work is proportional to
//!   the arithmetic and not to `n` (Gilbert–Peierls' result).
//! * **Threshold partial pivoting with diagonal preference**: the diagonal entry is the pivot unless it is
//!   smaller than `pivot_tol` times the largest candidate. For a structurally symmetric matrix and a symmetric
//!   ordering this keeps the factors symmetric in structure (so the fill is that of the ordering) while still
//!   protecting against tiny pivots. `pivot_tol = 1` is classical partial pivoting, `0` static diagonal pivoting.
//! * Generic over [`Scalar`] (real or complex) and independent of `faer`.
//!
//! Compared with `faer`'s supernodal LU this code is **not** BLAS-3 and has no supernodes: it wins when the
//! ordering matters more than the dense kernels (fill-reduction on 3-D graphs), and loses on matrices with
//! large dense fronts. `examples/lu_ordering_study.rs` and `docs/LEARNINGS.md` give the measurements.
//!
//! [`ShiftedSparseLu`] wraps it as a [`ShiftedSolver`], so the eigensolvers can use it wherever they use
//! [`ShiftedFactorization`](crate::linalg::direct::ShiftedFactorization). The ordering is computed once from the
//! union pattern; every shift repeats the (numeric + symbolic, as pivoting is dynamic) factorisation.

#![allow(clippy::needless_range_loop)]

use crate::linalg::csr::Csr;
use crate::linalg::direct::{FactorStats, ShiftedSolver, SolverError, Used, assemble};
use crate::linalg::linop::Preconditioner;
use crate::linalg::ordering::nested_dissection;
use crate::linalg::scalar::Scalar;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::SymbolicSparseColMatRef;
use faer::sparse::linalg::amd;
use std::time::Instant;

const NONE: usize = usize::MAX;

/// Column ordering `Q` of `P A Q = L U`.
#[derive(Debug, Clone, Default)]
pub enum LuOrdering {
    Natural,
    /// Approximate minimum degree (faer's implementation) on the pattern of `A + Aᵀ`.
    #[default]
    Amd,
    /// Level-set nested dissection of [`crate::linalg::ordering`] on the pattern of `A + Aᵀ`; leaf size in nodes.
    NestedDissection {
        leaf: usize,
    },
    /// `order[k]` = the original column eliminated `k`-th.
    Custom(Vec<usize>),
}

#[derive(Debug, Clone)]
pub struct LuOptions {
    pub ordering: LuOrdering,
    /// Diagonal pivot accepted if `|a_jj| ≥ pivot_tol · max_i |a_ij|` (over the not yet pivotal rows).
    pub pivot_tol: f64,
}

impl Default for LuOptions {
    fn default() -> Self {
        Self { ordering: LuOrdering::default(), pivot_tol: 0.1 }
    }
}

/// Sizes and cost of a factorisation.
#[derive(Debug, Clone, Copy, Default)]
pub struct LuStats {
    pub nnz_l: usize,
    pub nnz_u: usize,
    /// Scalar multiply–subtract operations of the numeric factorisation.
    pub flops: u64,
    pub order_secs: f64,
    pub factor_secs: f64,
    /// Pivots taken off the diagonal.
    pub off_diagonal_pivots: usize,
}

/// `P A Q = L U` with unit lower `L` (strictly lower part stored) and upper `U` (diagonal stored last in each column).
#[derive(Debug, Clone)]
pub struct SparseLu<T: Scalar> {
    n: usize,
    /// `q[k]` = original column in position `k`.
    q: Vec<usize>,
    /// `pinv[i]` = position of original row `i` (so `P b` is `y[pinv[i]] = b[i]`).
    pinv: Vec<usize>,
    l_ptr: Vec<usize>,
    l_idx: Vec<usize>,
    l_val: Vec<T>,
    u_ptr: Vec<usize>,
    u_idx: Vec<usize>,
    u_val: Vec<T>,
    stats: LuStats,
}

/// Pattern of `A + Aᵀ` (without the diagonal) from CSC arrays, as adjacency lists.
pub fn symmetric_adjacency(n: usize, col_ptr: &[u32], row_idx: &[u32]) -> Vec<Vec<u32>> {
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
    for j in 0..n {
        for &ri in &row_idx[col_ptr[j] as usize..col_ptr[j + 1] as usize] {
            let i = ri as usize;
            if i != j {
                adj[j].push(i as u32);
                adj[i].push(j as u32);
            }
        }
    }
    for l in &mut adj {
        l.sort_unstable();
        l.dedup();
    }
    adj
}

fn amd_order(n: usize, col_ptr: &[u32], row_idx: &[u32]) -> Result<Vec<usize>, SolverError> {
    if n == 0 {
        return Ok(Vec::new());
    }
    let pattern = SymbolicSparseColMatRef::new_checked(n, n, col_ptr, None, row_idx);
    let mut perm = vec![0u32; n];
    let mut perm_inv = vec![0u32; n];
    let mut mem = MemBuffer::new(amd::order_maybe_unsorted_scratch::<u32>(n, row_idx.len()));
    amd::order_maybe_unsorted(&mut perm, &mut perm_inv, pattern, amd::Control::default(), MemStack::new(&mut mem))
        .map_err(|_| SolverError::OutOfMemory)?;
    Ok(perm.into_iter().map(|p| p as usize).collect())
}

fn resolve_ordering(
    ordering: &LuOrdering,
    n: usize,
    col_ptr: &[u32],
    row_idx: &[u32],
) -> Result<Vec<usize>, SolverError> {
    let q = match ordering {
        LuOrdering::Natural => (0..n).collect(),
        LuOrdering::Amd => amd_order(n, col_ptr, row_idx)?,
        LuOrdering::NestedDissection { leaf } => {
            nested_dissection(&symmetric_adjacency(n, col_ptr, row_idx), (*leaf).max(2))
        }
        LuOrdering::Custom(o) => o.clone(),
    };
    let mut seen = vec![false; n];
    if q.len() != n || q.iter().any(|&c| c >= n || std::mem::replace(&mut seen[c], true)) {
        return Err(SolverError::InvalidArgument("ordering is not a permutation of 0..n"));
    }
    Ok(q)
}

impl<T: Scalar> SparseLu<T> {
    /// Factorises a square CSR matrix.
    pub fn factor(a: &Csr<T>, opts: &LuOptions) -> Result<Self, SolverError> {
        let n = a.nrows();
        if a.ncols() != n {
            return Err(SolverError::InvalidArgument("matrix must be square"));
        }
        let at = a.transpose(); // rows of Aᵀ = columns of A
        let col_ptr: Vec<u32> = at.indptr().iter().map(|&p| p as u32).collect();
        Self::factor_csc(n, &col_ptr, at.indices(), at.values(), opts)
    }

    /// Factorises a matrix given by its CSC arrays (`row_idx` need not be sorted).
    pub fn factor_csc(
        n: usize,
        col_ptr: &[u32],
        row_idx: &[u32],
        vals: &[T],
        opts: &LuOptions,
    ) -> Result<Self, SolverError> {
        let t0 = Instant::now();
        let q = resolve_ordering(&opts.ordering, n, col_ptr, row_idx)?;
        let order_secs = t0.elapsed().as_secs_f64();
        Self::factor_with_order(n, col_ptr, row_idx, vals, q, opts.pivot_tol, order_secs)
    }

    fn factor_with_order(
        n: usize,
        col_ptr: &[u32],
        row_idx: &[u32],
        vals: &[T],
        q: Vec<usize>,
        pivot_tol: f64,
        order_secs: f64,
    ) -> Result<Self, SolverError> {
        let t0 = Instant::now();
        let mut pinv = vec![NONE; n];
        let cap = col_ptr[n] as usize * 4 + n;
        let (mut l_ptr, mut u_ptr) = (Vec::with_capacity(n + 1), Vec::with_capacity(n + 1));
        let (mut l_idx, mut u_idx): (Vec<usize>, Vec<usize>) = (Vec::with_capacity(cap), Vec::with_capacity(cap));
        let (mut l_val, mut u_val): (Vec<T>, Vec<T>) = (Vec::with_capacity(cap), Vec::with_capacity(cap));
        l_ptr.push(0);
        u_ptr.push(0);
        let mut x = vec![T::ZERO; n];
        let mut mark = vec![0u32; n];
        let mut stamp = 0u32;
        let mut stack = vec![0usize; n];
        let mut pstack = vec![0usize; n];
        let mut xi = vec![0usize; n];
        let mut stats = LuStats { order_secs, ..Default::default() };

        for k in 0..n {
            let j = q[k];
            stamp += 1;
            // --- symbolic: reach of the pattern of A(:,j) in the graph of L ---------------------------------
            let mut top = n;
            for &r0 in &row_idx[col_ptr[j] as usize..col_ptr[j + 1] as usize] {
                let r0 = r0 as usize;
                if mark[r0] == stamp {
                    continue;
                }
                let mut head = 0isize;
                stack[0] = r0;
                while head >= 0 {
                    let h = head as usize;
                    let node = stack[h];
                    let jnew = pinv[node];
                    if mark[node] != stamp {
                        mark[node] = stamp;
                        pstack[h] = if jnew == NONE { 0 } else { l_ptr[jnew] };
                    }
                    let end = if jnew == NONE { 0 } else { l_ptr[jnew + 1] };
                    let mut finished = true;
                    for p in pstack[h]..end {
                        let child = l_idx[p];
                        if mark[child] == stamp {
                            continue;
                        }
                        pstack[h] = p + 1;
                        head += 1;
                        stack[head as usize] = child;
                        finished = false;
                        break;
                    }
                    if finished {
                        top -= 1;
                        xi[top] = node;
                        head -= 1;
                    }
                }
            }
            // --- numeric: x = L \ A(:,j) on the reach ---------------------------------------------------------
            for p in col_ptr[j] as usize..col_ptr[j + 1] as usize {
                x[row_idx[p] as usize] += vals[p];
            }
            for &i in &xi[top..n] {
                let r = pinv[i];
                if r == NONE {
                    continue;
                }
                let xi_val = x[i];
                if xi_val == T::ZERO {
                    continue;
                }
                for p in l_ptr[r]..l_ptr[r + 1] {
                    x[l_idx[p]] -= l_val[p] * xi_val;
                }
                stats.flops += (l_ptr[r + 1] - l_ptr[r]) as u64;
            }
            // --- pivot --------------------------------------------------------------------------------------------
            let mut best = NONE;
            let mut best_abs = -1.0f64;
            for &i in &xi[top..n] {
                if pinv[i] == NONE {
                    let a = x[i].modulus();
                    if a > best_abs {
                        best_abs = a;
                        best = i;
                    }
                }
            }
            if best == NONE || best_abs <= 0.0 || !best_abs.is_finite() {
                return Err(SolverError::Singular);
            }
            let mut piv = best;
            if pinv[j] == NONE && x[j].modulus() >= pivot_tol * best_abs && x[j].modulus() > 0.0 {
                piv = j;
            }
            if piv != j {
                stats.off_diagonal_pivots += 1;
            }
            let pivot = x[piv];
            for &i in &xi[top..n] {
                let r = pinv[i];
                if r != NONE {
                    u_idx.push(r);
                    u_val.push(x[i]);
                } else if i != piv {
                    l_idx.push(i); // original row index; mapped to positions at the end
                    l_val.push(x[i] / pivot);
                }
                x[i] = T::ZERO;
            }
            u_idx.push(k);
            u_val.push(pivot);
            pinv[piv] = k;
            l_ptr.push(l_idx.len());
            u_ptr.push(u_idx.len());
        }
        for i in &mut l_idx {
            *i = pinv[*i];
        }
        stats.nnz_l = l_idx.len();
        stats.nnz_u = u_idx.len();
        stats.factor_secs = t0.elapsed().as_secs_f64();
        Ok(Self { n, q, pinv, l_ptr, l_idx, l_val, u_ptr, u_idx, u_val, stats })
    }

    pub fn dim(&self) -> usize {
        self.n
    }

    pub fn stats(&self) -> &LuStats {
        &self.stats
    }

    /// Entries of `L + U` (the unit diagonal of `L` not counted).
    pub fn factor_nnz(&self) -> usize {
        self.stats.nnz_l + self.stats.nnz_u
    }

    /// The column ordering used.
    pub fn column_order(&self) -> &[usize] {
        &self.q
    }

    /// Solves `A x = b` in place.
    pub fn solve_in_place(&self, b: &mut [T]) {
        let n = self.n;
        assert_eq!(b.len(), n);
        let mut y = vec![T::ZERO; n];
        for (i, &bi) in b.iter().enumerate() {
            y[self.pinv[i]] = bi;
        }
        // L y = P b (unit lower, column oriented)
        for k in 0..n {
            let yk = y[k];
            if yk != T::ZERO {
                for p in self.l_ptr[k]..self.l_ptr[k + 1] {
                    y[self.l_idx[p]] -= self.l_val[p] * yk;
                }
            }
        }
        // U w = y (column oriented, diagonal stored last in each column)
        for k in (0..n).rev() {
            let end = self.u_ptr[k + 1];
            let wk = y[k] / self.u_val[end - 1];
            y[k] = wk;
            if wk != T::ZERO {
                for p in self.u_ptr[k]..end - 1 {
                    y[self.u_idx[p]] -= self.u_val[p] * wk;
                }
            }
        }
        for (k, &c) in self.q.iter().enumerate() {
            b[c] = y[k];
        }
    }

    /// Solves for `nrhs` right-hand sides stored as consecutive columns.
    pub fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) {
        for c in 0..nrhs {
            self.solve_in_place(&mut rhs[c * self.n..(c + 1) * self.n]);
        }
    }
}

impl<T: Scalar> Preconditioner<T> for SparseLu<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
        self.solve_in_place(z);
    }
}

/// [`ShiftedSolver`] on top of [`SparseLu`]: `(A − σB)⁻¹` with the ordering computed once.
pub struct ShiftedSparseLu<T: Scalar> {
    n: usize,
    col_ptr: Vec<u32>,
    row_idx: Vec<u32>,
    a: Vec<T>,
    b: Vec<T>,
    q: Vec<usize>,
    pivot_tol: f64,
    vals: Vec<T>,
    lu: Option<SparseLu<T>>,
    stats: FactorStats,
    solves: std::cell::Cell<usize>,
}

impl<T: Scalar> ShiftedSparseLu<T> {
    pub fn new(a: &Csr<T>, b: Option<&Csr<T>>, opts: &LuOptions) -> Result<Self, SolverError> {
        let t0 = Instant::now();
        let asm = assemble(a, b, false)?;
        let n = a.nrows();
        let q = resolve_ordering(&opts.ordering, n, &asm.col_ptr, &asm.row_idx)?;
        let nnz = asm.row_idx.len();
        Ok(Self {
            n,
            col_ptr: asm.col_ptr,
            row_idx: asm.row_idx,
            a: asm.a,
            b: asm.b,
            q,
            pivot_tol: opts.pivot_tol,
            vals: vec![T::ZERO; nnz],
            lu: None,
            stats: FactorStats {
                analyse_secs: t0.elapsed().as_secs_f64(),
                symbolic_count: 1,
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
            solves: std::cell::Cell::new(0),
        })
    }

    /// The factorisation of the current shift.
    pub fn factorization(&self) -> Option<&SparseLu<T>> {
        self.lu.as_ref()
    }
}

impl<T: Scalar> ShiftedSolver<T> for ShiftedSparseLu<T> {
    fn dim(&self) -> usize {
        self.n
    }

    fn set_shift(&mut self, sigma: T) -> Result<(), SolverError> {
        let t0 = Instant::now();
        self.vals.clear();
        self.vals.extend(self.a.iter().zip(&self.b).map(|(&x, &y)| x - sigma * y));
        self.stats.assemble_secs_last = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let lu = SparseLu::factor_with_order(
            self.n,
            &self.col_ptr,
            &self.row_idx,
            &self.vals,
            self.q.clone(),
            self.pivot_tol,
            0.0,
        )?;
        self.stats.numeric_secs_last = t1.elapsed().as_secs_f64();
        self.stats.numeric_secs_total += self.stats.numeric_secs_last;
        self.stats.numeric_count += 1;
        self.stats.used = Used::Lu;
        self.stats.factor_nnz = lu.factor_nnz();
        self.lu = Some(lu);
        Ok(())
    }

    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) {
        let lu = self.lu.as_ref().expect("set_shift must be called before solving");
        lu.solve_block_in_place(rhs, nrhs);
        self.solves.set(self.solves.get() + nrhs);
    }

    fn stats(&self) -> FactorStats {
        FactorStats { solves: self.solves.get(), ..self.stats }
    }
}
