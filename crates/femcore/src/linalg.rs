//! Linear algebra on top of faer: sparse direct solvers and small dense eigenproblems.
//!
//! `SparseSolver::new` factorises a real sparse matrix once, with a supernodal Cholesky factorisation
//! `A = L L^T` when it is symmetric positive definite and a sparse LU with partial pivoting otherwise, and then
//! solves for any number of right-hand sides (`solve_many` spreads independent solves over the thread pool). The
//! factorisation itself runs on the small factorisation pool of `par` and cannot be interrupted (the stage log
//! says so). `eigh` (symmetric) and `eig_general` (general, possibly complex eigenvalues) are dense routines
//! for the projected matrices of the Krylov and LOBPCG eigensolvers.

use crate::sparse::Csr;
use faer::linalg::solvers::Solve;
use faer::sparse::linalg::solvers::{Llt, Lu};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use faer::{Mat, Side};
use rayon::prelude::*;

/// Sparse direct solver: Cholesky for symmetric positive definite matrices, LU otherwise.
pub enum SparseSolver {
    Chol(Llt<usize, f64>, usize),
    Lu(Box<Lu<usize, f64>>, usize),
}

impl SparseSolver {
    /// Factorises `a`. With `spd` a Cholesky factorisation is tried first (only the lower
    /// triangle is read); if the matrix turns out not to be positive definite LU is used.
    pub fn new(a: &Csr, spd: bool) -> Result<Self, String> {
        assert_eq!(a.nrows, a.ncols);
        let n = a.nrows;
        if n == 0 {
            return Err("empty matrix".into());
        }
        crate::progress::check()?;
        let _stage = crate::progress::stage(format!("factorising a sparse {} matrix ({n} unknowns, {} non-zeros); this cannot be interrupted", if spd { "symmetric" } else { "general" }, a.nnz()));
        // CSR of A^T is CSC of A
        let t = a.transpose();
        let sym = SymbolicSparseColMatRef::new_checked(n, n, &t.indptr, None, &t.indices);
        let mat = SparseColMatRef::new(sym, &t.data);
        crate::par::on_factor_pool(|| {
            if spd {
                if let Ok(f) = mat.sp_cholesky(Side::Lower) {
                    return Ok(SparseSolver::Chol(f, n));
                }
            }
            match mat.sp_lu() {
                Ok(f) => Ok(SparseSolver::Lu(Box::new(f), n)),
                Err(e) => Err(format!("sparse factorisation failed (singular matrix?): {e:?}")),
            }
        })
    }

    pub fn n(&self) -> usize {
        match self {
            SparseSolver::Chol(_, n) | SparseSolver::Lu(_, n) => *n,
        }
    }

    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n();
        assert_eq!(b.len(), n);
        let mut rhs = Mat::<f64>::from_fn(n, 1, |i, _| b[i]);
        match self {
            SparseSolver::Chol(f, _) => f.solve_in_place(rhs.as_mut()),
            SparseSolver::Lu(f, _) => f.solve_in_place(rhs.as_mut()),
        }
        (0..n).map(|i| rhs[(i, 0)]).collect()
    }

    /// Several right-hand sides, the independent triangular solves spread over threads.
    pub fn solve_many(&self, bs: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = self.n();
        let nt = if bs.len() > 1 && n > 5_000 { crate::par::threads().min(bs.len()) } else { 1 };
        if nt <= 1 {
            return bs.iter().map(|b| self.solve(b)).collect();
        }
        let chunk = bs.len().div_ceil(nt);
        let parts: Vec<Vec<Vec<f64>>> = bs.par_chunks(chunk).map(|cols| cols.iter().map(|b| self.solve(b)).collect::<Vec<Vec<f64>>>()).collect();
        parts.into_iter().flatten().collect()
    }
}

/// Dense matrix as row-major `Vec<f64>` with `n` columns.
pub fn to_faer(a: &[f64], nrows: usize, ncols: usize) -> Mat<f64> {
    Mat::from_fn(nrows, ncols, |i, j| a[i * ncols + j])
}

/// Eigen decomposition of a symmetric dense matrix (row-major, only the lower triangle is
/// read). Returns the ascending eigenvalues and the eigenvectors as columns
/// (`vecs[j]` is the j-th eigenvector).
pub fn eigh(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let m = to_faer(a, n, n);
    let e = m.self_adjoint_eigen(Side::Lower).map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    let u = e.U();
    let s = e.S();
    let vals: Vec<f64> = (0..n).map(|i| s.column_vector()[i]).collect();
    let vecs: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| u[(i, j)]).collect()).collect();
    Ok((vals, vecs))
}

/// Generalised symmetric-definite problem `K x = lambda M x` (dense, row-major, M positive
/// definite) via the Cholesky factor of M. Returns ascending eigenvalues and M-orthonormal
/// eigenvectors.
pub fn gen_eigh(k: &[f64], m: &[f64], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let mm = to_faer(m, n, n);
    let llt = mm.llt(Side::Lower).map_err(|e| format!("mass matrix is not positive definite: {e:?}"))?;
    let l = llt.L();
    // C = L^-1 K L^-T
    let mut c = to_faer(k, n, n);
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, c.as_mut(), faer::get_global_parallelism());
    let mut ct = c.transpose().to_owned();
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, ct.as_mut(), faer::get_global_parallelism());
    // symmetrise against round-off
    let cs = Mat::from_fn(n, n, |i, j| 0.5 * (ct[(i, j)] + ct[(j, i)]));
    let e = cs.self_adjoint_eigen(Side::Lower).map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    let (u, s) = (e.U(), e.S());
    let vals: Vec<f64> = (0..n).map(|i| s.column_vector()[i]).collect();
    // x = L^-T y
    let mut y = Mat::from_fn(n, n, |i, j| u[(i, j)]);
    faer::linalg::triangular_solve::solve_upper_triangular_in_place(
        l.transpose(),
        y.as_mut(),
        faer::get_global_parallelism(),
    );
    let vecs: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| y[(i, j)]).collect()).collect();
    Ok((vals, vecs))
}

/// One eigenpair of a general real matrix.
pub struct GeneralEig {
    pub value: (f64, f64),
    pub vector: Vec<(f64, f64)>,
}

/// Eigenpairs of a general real matrix (row-major).
pub fn eig_general(a: &[f64], n: usize) -> Result<Vec<GeneralEig>, String> {
    let m = to_faer(a, n, n);
    let e = m.eigen().map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    let (u, s) = (e.U(), e.S());
    Ok((0..n)
        .map(|j| {
            let z = s.column_vector()[j];
            GeneralEig {
                value: (z.re, z.im),
                vector: (0..n).map(|i| (u[(i, j)].re, u[(i, j)].im)).collect(),
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Coo;

    fn laplace_1d(n: usize) -> Csr {
        let mut c = Coo::new(n, n);
        for i in 0..n {
            c.push(i, i, 2.0);
            if i > 0 {
                c.push(i, i - 1, -1.0);
                c.push(i - 1, i, -1.0);
            }
        }
        c.to_csr()
    }

    #[test]
    fn sparse_solve_chol_and_lu() {
        let a = laplace_1d(50);
        let b: Vec<f64> = (0..50).map(|i| (i as f64).sin()).collect();
        for spd in [true, false] {
            let s = SparseSolver::new(&a, spd).unwrap();
            let x = s.solve(&b);
            let r = a.matvec(&x);
            let err = r.iter().zip(&b).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
            assert!(err < 1e-12, "{err}");
        }
    }

    #[test]
    fn dense_eigs() {
        let n = 20;
        let a = laplace_1d(n).to_dense();
        let (v, _) = eigh(&a, n).unwrap();
        for (k, &lam) in v.iter().enumerate() {
            let th = std::f64::consts::PI * (k + 1) as f64 / (n + 1) as f64;
            assert!((lam - (2.0 - 2.0 * th.cos())).abs() < 1e-12);
        }
        let mut m = vec![0.0; n * n];
        for i in 0..n {
            m[i * n + i] = 2.0;
        }
        let (v2, vecs) = gen_eigh(&a, &m, n).unwrap();
        for k in 0..n {
            assert!((v2[k] - 0.5 * v[k]).abs() < 1e-12);
            // M-orthonormal
            let nrm: f64 = vecs[k].iter().map(|x| 2.0 * x * x).sum();
            assert!((nrm - 1.0).abs() < 1e-10);
        }
        // general: a rotation-like matrix has eigenvalues +-i
        let g = eig_general(&[0.0, -1.0, 1.0, 0.0], 2).unwrap();
        assert!(g.iter().all(|e| e.value.0.abs() < 1e-14 && (e.value.1.abs() - 1.0).abs() < 1e-14));
    }
}
