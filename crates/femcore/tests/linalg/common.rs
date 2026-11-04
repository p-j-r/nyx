#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]
//! Shared helpers for the integration tests: independent residual checks and analytic spectra.

use femcore::linalg::{C64, Csr, LinearOperator, Scalar};

pub fn norm<T: Scalar>(v: &[T]) -> f64 {
    v.iter().map(|x| x.modulus_sq()).sum::<f64>().sqrt()
}

/// `‖A x − λ B x‖ / ‖x‖` (B = identity when `None`), computed with plain loops over the
/// matrices so it shares no code with the solvers.
pub fn true_residual<T: Scalar>(a: &Csr<T>, b: Option<&Csr<T>>, lambda: C64, x: &[T]) -> f64 {
    let n = a.nrows();
    let mut ax = vec![T::ZERO; n];
    a.spmv(x, &mut ax);
    let mut bx = x.to_vec();
    if let Some(b) = b {
        b.spmv(x, &mut bx);
    }
    let mut s = 0.0;
    for i in 0..n {
        let d = ax[i].to_c64() - lambda * bx[i].to_c64();
        s += d.norm_sqr();
    }
    s.sqrt() / norm(x)
}

/// Eigenvalues of `tridiag(-1, 2, -1)` of order `n`, ascending.
pub fn laplacian_1d_eigenvalues(n: usize) -> Vec<f64> {
    (1..=n).map(|k| 2.0 - 2.0 * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos()).collect()
}

/// Eigenvalues of the 5-point Laplacian on an `nx × ny` grid, ascending (with multiplicity).
pub fn laplacian_2d_eigenvalues(nx: usize, ny: usize) -> Vec<f64> {
    let (ex, ey) = (laplacian_1d_eigenvalues(nx), laplacian_1d_eigenvalues(ny));
    let mut v: Vec<f64> = ex.iter().flat_map(|a| ey.iter().map(move |b| a + b)).collect();
    v.sort_by(f64::total_cmp);
    v
}

/// Dense operator from a row-major matrix (tests with tiny problems).
pub struct DenseOp<T: Scalar> {
    pub n: usize,
    pub a: Vec<T>,
}

impl<T: Scalar> LinearOperator<T> for DenseOp<T> {
    fn dim(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        for i in 0..self.n {
            y[i] = (0..self.n).fold(T::ZERO, |s, j| s + self.a[i * self.n + j] * x[j]);
        }
    }
}

/// Solves `B X = A` densely (Gaussian elimination with partial pivoting) and returns the
/// row-major `B⁻¹ A`.
pub fn dense_b_inv_a(a: &Csr<f64>, b: &Csr<f64>) -> Vec<f64> {
    let n = a.nrows();
    let (mut m, mut x) = (b.to_dense(), a.to_dense());
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| m[i * n + k].abs().total_cmp(&m[j * n + k].abs())).unwrap();
        if p != k {
            for c in 0..n {
                m.swap(k * n + c, p * n + c);
                x.swap(k * n + c, p * n + c);
            }
        }
        let piv = m[k * n + k];
        for i in k + 1..n {
            let f = m[i * n + k] / piv;
            for c in 0..n {
                let (mk, xk) = (m[k * n + c], x[k * n + c]);
                m[i * n + c] -= f * mk;
                x[i * n + c] -= f * xk;
            }
        }
    }
    for c in 0..n {
        for i in (0..n).rev() {
            let mut s = x[i * n + c];
            for j in i + 1..n {
                s -= m[i * n + j] * x[j * n + c];
            }
            x[i * n + c] = s / m[i * n + i];
        }
    }
    x
}
