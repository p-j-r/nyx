//! Small dense linear algebra for the projected (Rayleigh–Ritz / Schur) problems of the
//! Krylov methods, the role of SLEPc's `DS`.
//!
//! * [`hermitian_eig`]: Hermitian eigendecomposition through `faer` (see `faer_glue.rs`).
//! * [`schur`]: complex Schur form by Hessenberg reduction and shifted QR, written here because
//!   `faer` does not expose Schur reordering and the Krylov–Schur restart needs it.
//! * [`reorder`]: adjacent-swap reordering of the Schur form, [`triangular_eigenvector`].

mod schur;

pub use schur::{reorder, schur, triangular_eigenvector};

use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::scalar::Scalar;
use std::ops::{Index, IndexMut};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DenseError {
    /// The QR iteration or the dense eigensolver did not converge.
    NoConvergence,
    NonFinite,
}

impl std::fmt::Display for DenseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DenseError::NoConvergence => write!(f, "dense eigen/Schur iteration did not converge"),
            DenseError::NonFinite => write!(f, "dense problem contains NaN or infinity"),
        }
    }
}

impl std::error::Error for DenseError {}

/// Column-major dense matrix with leading dimension equal to the row count.
#[derive(Debug, Clone, PartialEq)]
pub struct Dense<T> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
}

impl<T: Scalar> Dense<T> {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self { rows, cols, data: vec![T::ZERO; rows * cols] }
    }
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = T::ONE;
        }
        m
    }
    pub fn from_fn(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> T) -> Self {
        let mut m = Self::zeros(rows, cols);
        for j in 0..cols {
            for i in 0..rows {
                m[(i, j)] = f(i, j);
            }
        }
        m
    }
    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }
    pub fn col(&self, j: usize) -> &[T] {
        &self.data[j * self.rows..(j + 1) * self.rows]
    }
    pub fn col_mut(&mut self, j: usize) -> &mut [T] {
        &mut self.data[j * self.rows..(j + 1) * self.rows]
    }
    pub fn is_finite(&self) -> bool {
        self.data.iter().all(|v| v.is_finite())
    }
    pub fn frobenius(&self) -> f64 {
        self.data.iter().map(|v| v.modulus_sq()).sum::<f64>().sqrt()
    }
    pub fn adjoint(&self) -> Self {
        Self::from_fn(self.cols, self.rows, |i, j| self[(j, i)].conj())
    }
    /// Matrix product `self · rhs`.
    pub fn matmul(&self, rhs: &Dense<T>) -> Dense<T> {
        assert_eq!(self.cols, rhs.rows);
        let mut out = Dense::zeros(self.rows, rhs.cols);
        for j in 0..rhs.cols {
            for k in 0..self.cols {
                let b = rhs[(k, j)];
                if b != T::ZERO {
                    for i in 0..self.rows {
                        let v = self[(i, k)] * b;
                        out[(i, j)] += v;
                    }
                }
            }
        }
        out
    }
}

impl<T> Index<(usize, usize)> for Dense<T> {
    type Output = T;
    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &T {
        debug_assert!(i < self.rows && j < self.cols);
        &self.data[j * self.rows + i]
    }
}

impl<T> IndexMut<(usize, usize)> for Dense<T> {
    #[inline]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        debug_assert!(i < self.rows && j < self.cols);
        &mut self.data[j * self.rows + i]
    }
}

/// Eigendecomposition `A = U diag(w) Uᴴ` of a Hermitian matrix (only the lower triangle is
/// read). Eigenvalues are returned in ascending order with the matching orthonormal columns.
pub fn hermitian_eig<T: FaerScalar>(a: &Dense<T>) -> Result<(Vec<f64>, Dense<T>), DenseError> {
    crate::linalg::faer_glue::hermitian_eig(a)
}

/// Solves the small dense complex system `A x = b` by Gaussian elimination with partial pivoting
/// (used for tiny projected systems such as the deflation Gram matrices; `A` is `n × n`).
/// Returns `None` if `A` is numerically singular.
pub fn solve_small(a: &Dense<crate::linalg::scalar::C64>, b: &[crate::linalg::scalar::C64]) -> Option<Vec<crate::linalg::scalar::C64>> {
    let n = a.rows();
    assert!(a.cols() == n && b.len() == n);
    let mut m = a.clone();
    let mut x = b.to_vec();
    let scale = m.frobenius().max(f64::MIN_POSITIVE);
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| m[(i, k)].norm().total_cmp(&m[(j, k)].norm()))?;
        if m[(p, k)].norm() <= 1e-14 * scale {
            return None;
        }
        if p != k {
            for c in 0..n {
                let t = m[(k, c)];
                m[(k, c)] = m[(p, c)];
                m[(p, c)] = t;
            }
            x.swap(k, p);
        }
        for i in k + 1..n {
            let f = m[(i, k)] / m[(k, k)];
            for c in k..n {
                let v = m[(k, c)];
                m[(i, c)] -= f * v;
            }
            let xk = x[k];
            x[i] -= f * xk;
        }
    }
    for i in (0..n).rev() {
        let mut sum = x[i];
        for c in i + 1..n {
            sum -= m[(i, c)] * x[c];
        }
        x[i] = sum / m[(i, i)];
    }
    Some(x)
}
