//! **Mixed-precision shifted factorisation** with iterative refinement (port item S3): the sparse LU of `A − σB`
//! is computed in single precision (`f32` for real, `Complex<f32>` for complex matrices) — about half the
//! factor memory and a faster factorisation — and every solve is refined in double precision,
//! `x ← x + LU⁻¹(b − (A − σB)x)`, until the residual is at the double-precision level.
//!
//! This only works while `κ(A − σB)·ε₃₂ < 1`. A shifted matrix close to a singularity (exactly what shift-and-invert
//! eigensolvers produce near a converged eigenvalue) breaks that: the refinement then stagnates or diverges. The solver
//! **detects it** (the residual must drop by at least a factor 2 per step) and falls back to a double-precision
//! factorisation of the same matrix, built on first need, so results stay correct and only the speed advantage is lost.
//! [`MixedStats`] reports how often refinement was used and how many steps it took.
//!
//! Implements [`ShiftedSolver`], so it can replace [`ShiftedFactorization`] in a [`ShiftInvertOperator`](crate::solvers::direct::ShiftInvertOperator).

use crate::solvers::csr::Csr;
use crate::solvers::direct::{
    Definiteness, FactorStats, ShiftedFactorization, ShiftedSolver, SolverError, ThreadPolicy, Used, assemble,
};
use crate::solvers::faer_glue::{FaerScalar, clear_upper_vector_state};
use crate::solvers::scalar::C64;
use crate::solvers::vecops::norm2;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::SparseColMatRef;
use faer::sparse::SymbolicSparseColMatRef;
use faer::sparse::linalg::lu::{LuRef, NumericLu, SymbolicLu, factorize_symbolic_lu};
use faer::traits::ComplexField;
use faer::{Conj, MatMut};
use std::cell::{Cell, RefCell};
use std::time::Instant;

/// Scalars with a single-precision counterpart usable by `faer`.
pub trait LowPrecision: FaerScalar {
    type Low: ComplexField + Copy;
    fn to_low(self) -> Self::Low;
    fn from_low(l: Self::Low) -> Self;
}

impl LowPrecision for f64 {
    type Low = f32;
    fn to_low(self) -> f32 {
        self as f32
    }
    fn from_low(l: f32) -> f64 {
        l as f64
    }
}

impl LowPrecision for C64 {
    type Low = faer::c32;
    fn to_low(self) -> faer::c32 {
        faer::c32::new(self.re as f32, self.im as f32)
    }
    fn from_low(l: faer::c32) -> C64 {
        C64::new(l.re as f64, l.im as f64)
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct MixedStats {
    /// Right-hand sides solved through refinement / by the double-precision fallback.
    pub refined_solves: usize,
    pub fallback_solves: usize,
    /// Total refinement steps (one single-precision solve each).
    pub refinement_steps: usize,
    /// Single-precision factorisations that failed (singular in `f32`) and numeric factorisations done.
    pub low_factor_failures: usize,
    pub low_factor_secs: f64,
}

pub struct MixedPrecisionLu<T: LowPrecision> {
    n: usize,
    a: Csr<T>,
    b: Option<Csr<T>>,
    col_ptr: Vec<u32>,
    row_idx: Vec<u32>,
    a_vals: Vec<T>,
    b_vals: Vec<T>,
    sym: Option<SymbolicLu<u32>>,
    num: NumericLu<u32, T::Low>,
    low_ok: bool,
    sigma: T,
    /// Required drop of the residual per refinement step.
    pub min_contraction: f64,
    pub max_refine: usize,
    /// Relative residual to reach (default 1e-13).
    pub tol: f64,
    /// Parallelism of the single-precision factorisation and solves (same defaults as [`ShiftedFactorization`]).
    pub threads: ThreadPolicy,
    exact: RefCell<Option<ShiftedFactorization<T>>>,
    stats: Cell<MixedStats>,
    base: FactorStats,
}

impl<T: LowPrecision> MixedPrecisionLu<T> {
    pub fn new(a: &Csr<T>, b: Option<&Csr<T>>) -> Result<Self, SolverError> {
        let asm = assemble(a, b, false)?;
        Ok(Self {
            n: a.nrows(),
            a: a.clone(),
            b: b.cloned(),
            col_ptr: asm.col_ptr,
            row_idx: asm.row_idx,
            a_vals: asm.a,
            b_vals: asm.b,
            sym: None,
            num: NumericLu::new(),
            low_ok: false,
            sigma: T::ZERO,
            min_contraction: 0.5,
            max_refine: 12,
            tol: 1e-13,
            threads: ThreadPolicy::default(),
            exact: RefCell::new(None),
            stats: Cell::new(MixedStats::default()),
            base: FactorStats {
                analyse_secs: 0.0,
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

    pub fn mixed_stats(&self) -> MixedStats {
        self.stats.get()
    }

    fn pattern(&self) -> SymbolicSparseColMatRef<'_, u32> {
        SymbolicSparseColMatRef::new_checked(self.n, self.n, &self.col_ptr, None, &self.row_idx)
    }

    /// `y ← (A − σB) x` in working precision.
    fn apply_shifted(&self, x: &[T], y: &mut [T]) {
        self.a.spmv(x, y);
        if let Some(b) = &self.b {
            let mut t = vec![T::ZERO; self.n];
            b.spmv(x, &mut t);
            for (yi, ti) in y.iter_mut().zip(&t) {
                *yi -= self.sigma * *ti;
            }
        } else {
            for (yi, xi) in y.iter_mut().zip(x) {
                *yi -= self.sigma * *xi;
            }
        }
    }

    fn solve_low(&self, v: &mut [T::Low]) {
        let sym = self.sym.as_ref().unwrap();
        let par = self.threads.solve.par();
        let mut mem = MemBuffer::new(sym.solve_in_place_scratch::<T::Low>(1, par));
        let m = MatMut::from_column_major_slice_mut(v, self.n, 1);
        LuRef::new_unchecked(sym, &self.num).solve_in_place_with_conj(Conj::No, m, par, MemStack::new(&mut mem));
        clear_upper_vector_state();
    }

    fn exact_solve(&self, col: &mut [T]) {
        let mut slot = self.exact.borrow_mut();
        if slot.is_none() {
            let mut f = ShiftedFactorization::new(&self.a, self.b.as_ref(), Definiteness::General).expect("pattern");
            f.set_shift(self.sigma).expect("double-precision factorisation");
            *slot = Some(f);
        }
        slot.as_ref().unwrap().solve_block_in_place(col, 1);
        let mut st = self.stats.get();
        st.fallback_solves += 1;
        self.stats.set(st);
    }

    /// One right-hand side by refinement; `false` if the refinement did not converge.
    fn refine(&self, rhs: &[T], x: &mut [T]) -> bool {
        let n = self.n;
        let bn = norm2(rhs);
        x.fill(T::ZERO);
        if bn == 0.0 {
            return true;
        }
        let mut r = rhs.to_vec();
        let mut prev = bn;
        let mut ax = vec![T::ZERO; n];
        let mut steps = 0;
        let mut ok = false;
        for _ in 0..self.max_refine {
            let rn = norm2(&r);
            let scale = 1.0 / rn;
            let mut low: Vec<T::Low> = r.iter().map(|v| v.scale(scale).to_low()).collect();
            self.solve_low(&mut low);
            for (xi, l) in x.iter_mut().zip(&low) {
                *xi += T::from_low(*l).scale(rn);
            }
            steps += 1;
            self.apply_shifted(x, &mut ax);
            for i in 0..n {
                r[i] = rhs[i] - ax[i];
            }
            let new = norm2(&r);
            if !new.is_finite() {
                break;
            }
            if new <= self.tol * bn {
                ok = true;
                break;
            }
            if new > self.min_contraction * prev {
                break; // stagnation or divergence: single precision is not enough for this matrix
            }
            prev = new;
        }
        let mut st = self.stats.get();
        st.refinement_steps += steps;
        if ok {
            st.refined_solves += 1;
        }
        self.stats.set(st);
        ok
    }
}

impl<T: LowPrecision> ShiftedSolver<T> for MixedPrecisionLu<T> {
    fn dim(&self) -> usize {
        self.n
    }

    fn set_shift(&mut self, sigma: T) -> Result<(), SolverError> {
        self.sigma = sigma;
        *self.exact.borrow_mut() = None;
        if self.sym.is_none() {
            let t0 = Instant::now();
            self.sym =
                Some(factorize_symbolic_lu(self.pattern(), Default::default()).map_err(|_| SolverError::OutOfMemory)?);
            self.base.analyse_secs += t0.elapsed().as_secs_f64();
            self.base.symbolic_count += 1;
        }
        let vals: Vec<T::Low> = self.a_vals.iter().zip(&self.b_vals).map(|(&a, &b)| (a - sigma * b).to_low()).collect();
        let t0 = Instant::now();
        let sym = self.sym.as_ref().unwrap();
        let par = self.threads.factor.par();
        let mut mem = MemBuffer::new(sym.factorize_numeric_lu_scratch::<T::Low>(par, Default::default()));
        // (the pattern is built from the fields directly so that `self.num` can be borrowed mutably)
        let pattern = SymbolicSparseColMatRef::new_checked(self.n, self.n, &self.col_ptr, None, &self.row_idx);
        let mat = SparseColMatRef::new(pattern, &vals);
        let res =
            sym.factorize_numeric_lu::<T::Low>(&mut self.num, mat, par, MemStack::new(&mut mem), Default::default());
        clear_upper_vector_state();
        let secs = t0.elapsed().as_secs_f64();
        self.low_ok = res.is_ok();
        let mut st = self.stats.get();
        st.low_factor_secs += secs;
        if !self.low_ok {
            st.low_factor_failures += 1;
        }
        self.stats.set(st);
        self.base.numeric_secs_last = secs;
        self.base.numeric_secs_total += secs;
        self.base.numeric_count += 1;
        self.base.used = Used::Lu;
        Ok(())
    }

    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) {
        let n = self.n;
        assert!(rhs.len() >= n * nrhs);
        for c in 0..nrhs {
            let col = &mut rhs[c * n..(c + 1) * n];
            if self.low_ok {
                let b = col.to_vec();
                if self.refine(&b, col) {
                    continue;
                }
                col.copy_from_slice(&b);
            }
            self.exact_solve(col);
        }
    }

    fn stats(&self) -> FactorStats {
        self.base
    }
}
