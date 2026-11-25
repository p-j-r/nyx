//! Compressed sparse row matrix over a generic [`Scalar`], plus test-problem generators.
//!
//! nyx has its own assembly types (`CCoo`/CSC); this type is the minimal interchange format
//! used by the solvers, the Matrix Market export and the tests. Convert with
//! [`Csr::from_triplets`] or [`Csr::try_new`].

use crate::linalg::linop::LinearOperator;
use crate::linalg::par;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::CHUNK;
use std::fmt;

mod assembly;
pub use assembly::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CsrError {
    BadIndptr(&'static str),
    LengthMismatch,
    ColumnOutOfRange { row: usize, col: usize },
    DimensionTooLarge,
}

impl fmt::Display for CsrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CsrError::BadIndptr(m) => write!(f, "invalid indptr: {m}"),
            CsrError::LengthMismatch => write!(f, "indices/values length does not match indptr"),
            CsrError::ColumnOutOfRange { row, col } => write!(f, "column {col} out of range in row {row}"),
            CsrError::DimensionTooLarge => write!(f, "dimension exceeds u32 index range"),
        }
    }
}

impl std::error::Error for CsrError {}

/// CSR matrix with `u32` column indices.
///
/// The invariants checked in [`try_new`](Self::try_new) (monotone `indptr`, every column index
/// `< ncols`) justify the unchecked indexing in the SpMV kernel; the fields are private so they
/// cannot be broken afterwards.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr<T: Scalar> {
    nrows: usize,
    ncols: usize,
    indptr: Vec<usize>,
    indices: Vec<u32>,
    values: Vec<T>,
}

impl<T: Scalar> Csr<T> {
    pub fn try_new(
        nrows: usize,
        ncols: usize,
        indptr: Vec<usize>,
        indices: Vec<u32>,
        values: Vec<T>,
    ) -> Result<Self, CsrError> {
        if ncols > u32::MAX as usize {
            return Err(CsrError::DimensionTooLarge);
        }
        let Some(expected_indptr_len) = nrows.checked_add(1) else {
            return Err(CsrError::DimensionTooLarge);
        };
        if indptr.len() != expected_indptr_len {
            return Err(CsrError::BadIndptr("length must be nrows + 1"));
        }
        if indptr[0] != 0 {
            return Err(CsrError::BadIndptr("must start at 0"));
        }
        if indptr.windows(2).any(|w| w[0] > w[1]) {
            return Err(CsrError::BadIndptr("must be non-decreasing"));
        }
        let nnz = indptr[nrows];
        if indices.len() != nnz || values.len() != nnz {
            return Err(CsrError::LengthMismatch);
        }
        for row in 0..nrows {
            for &c in &indices[indptr[row]..indptr[row + 1]] {
                if c as usize >= ncols {
                    return Err(CsrError::ColumnOutOfRange { row, col: c as usize });
                }
            }
        }
        Ok(Self { nrows, ncols, indptr, indices, values })
    }

    /// Builds from `(row, col, value)` triplets; duplicates are summed, columns sorted.
    pub fn from_triplets(nrows: usize, ncols: usize, triplets: &[(usize, usize, T)]) -> Result<Self, CsrError> {
        if nrows.checked_add(1).is_none() || ncols > u32::MAX as usize {
            return Err(CsrError::DimensionTooLarge);
        }
        let mut count = vec![0usize; nrows + 1];
        for &(r, c, _) in triplets {
            if r >= nrows {
                return Err(CsrError::BadIndptr("row index out of range"));
            }
            if c >= ncols {
                return Err(CsrError::ColumnOutOfRange { row: r, col: c });
            }
            count[r + 1] += 1;
        }
        for i in 0..nrows {
            count[i + 1] += count[i];
        }
        let mut next = count.clone();
        let mut entries = vec![(0u32, T::ZERO); triplets.len()];
        for &(r, c, v) in triplets {
            let cu = u32::try_from(c).map_err(|_| CsrError::DimensionTooLarge)?;
            entries[next[r]] = (cu, v);
            next[r] += 1;
        }
        let mut indptr = Vec::with_capacity(nrows + 1);
        let mut indices = Vec::with_capacity(triplets.len());
        let mut values = Vec::with_capacity(triplets.len());
        indptr.push(0);
        for r in 0..nrows {
            let row = &mut entries[count[r]..count[r + 1]];
            row.sort_by_key(|&(c, _)| c);
            let mut last: Option<u32> = None;
            for &(c, v) in row.iter() {
                if last == Some(c) {
                    *values.last_mut().unwrap() += v;
                } else {
                    indices.push(c);
                    values.push(v);
                    last = Some(c);
                }
            }
            indptr.push(indices.len());
        }
        Self::try_new(nrows, ncols, indptr, indices, values)
    }

    pub fn identity(n: usize) -> Self {
        assert!(n <= u32::MAX as usize, "Csr::identity dimension {n} overflows u32");
        let indices: Vec<u32> = (0..n).map(|i| u32::try_from(i).expect("u32 index")).collect();
        Self::try_new(n, n, (0..=n).collect(), indices, vec![T::ONE; n]).unwrap()
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }
    pub fn ncols(&self) -> usize {
        self.ncols
    }
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
    pub fn indptr(&self) -> &[usize] {
        &self.indptr
    }
    pub fn indices(&self) -> &[u32] {
        &self.indices
    }
    pub fn values(&self) -> &[T] {
        &self.values
    }

    /// `(columns, values)` of one row.
    pub fn row(&self, i: usize) -> (&[u32], &[T]) {
        let r = self.indptr[i]..self.indptr[i + 1];
        (&self.indices[r.clone()], &self.values[r])
    }

    /// Columns are strictly increasing within every row.
    pub fn has_sorted_unique_rows(&self) -> bool {
        (0..self.nrows).all(|i| self.row(i).0.windows(2).all(|w| w[0] < w[1]))
    }

    /// Diagonal entries (zero where structurally absent).
    pub fn diagonal(&self) -> Vec<T> {
        (0..self.nrows.min(self.ncols))
            .map(|i| {
                let (c, v) = self.row(i);
                c.iter().zip(v).filter(|&(&c, _)| c as usize == i).fold(T::ZERO, |a, (_, &v)| a + v)
            })
            .collect()
    }

    /// Plain transpose `Aᵀ` (no conjugation).
    pub fn transpose(&self) -> Csr<T> {
        assert!(self.nrows <= u32::MAX as usize, "Csr::transpose nrows overflows u32");
        let mut count = vec![0usize; self.ncols + 1];
        for &c in &self.indices {
            count[c as usize + 1] += 1;
        }
        for i in 0..self.ncols {
            count[i + 1] += count[i];
        }
        let mut next = count.clone();
        let mut indices = vec![0u32; self.nnz()];
        let mut values = vec![T::ZERO; self.nnz()];
        for i in 0..self.nrows {
            let iu = u32::try_from(i).expect("u32 row index");
            let (cols, vals) = self.row(i);
            for (&c, &v) in cols.iter().zip(vals) {
                let p = next[c as usize];
                next[c as usize] += 1;
                indices[p] = iu;
                values[p] = v;
            }
        }
        Self::try_new(self.ncols, self.nrows, count, indices, values).expect("valid transposed CSR")
    }

    /// `A = Aᵀ` (complex symmetric for complex scalars) up to `tol·max|a|`, pattern included.
    pub fn is_symmetric(&self, tol: f64) -> bool {
        self.matches_transpose(tol, false)
    }

    /// `A = Aᴴ` up to `tol·max|a|`.
    pub fn is_hermitian(&self, tol: f64) -> bool {
        self.matches_transpose(tol, true)
    }

    fn matches_transpose(&self, tol: f64, conj: bool) -> bool {
        if self.nrows != self.ncols {
            return false;
        }
        let t = self.transpose();
        if t.indptr != self.indptr || t.indices != self.indices {
            return false;
        }
        let scale = self.values.iter().fold(0.0f64, |m, v| m.max(v.modulus()));
        self.values.iter().zip(&t.values).all(|(a, b)| (*a - if conj { b.conj() } else { *b }).modulus() <= tol * scale)
    }

    /// Dense row-major copy (tests and tiny problems only).
    pub fn to_dense(&self) -> Vec<T> {
        let mut d = vec![T::ZERO; self.nrows * self.ncols];
        for i in 0..self.nrows {
            let (c, v) = self.row(i);
            for (&c, &v) in c.iter().zip(v) {
                d[i * self.ncols + c as usize] += v;
            }
        }
        d
    }

    /// `α·A + β·B` on the union of the two patterns.
    pub fn axpby(alpha: T, a: &Csr<T>, beta: T, b: &Csr<T>) -> Csr<T> {
        assert!(a.nrows == b.nrows && a.ncols == b.ncols);
        let mut indptr = Vec::with_capacity(a.nrows + 1);
        let (mut indices, mut values) = (Vec::new(), Vec::new());
        indptr.push(0);
        for i in 0..a.nrows {
            let ((ca, va), (cb, vb)) = (a.row(i), b.row(i));
            let (mut p, mut q) = (0, 0);
            while p < ca.len() || q < cb.len() {
                let (col, val);
                if q >= cb.len() || (p < ca.len() && ca[p] < cb[q]) {
                    (col, val) = (ca[p], alpha * va[p]);
                    p += 1;
                } else if p >= ca.len() || cb[q] < ca[p] {
                    (col, val) = (cb[q], beta * vb[q]);
                    q += 1;
                } else {
                    (col, val) = (ca[p], alpha * va[p] + beta * vb[q]);
                    p += 1;
                    q += 1;
                }
                indices.push(col);
                values.push(val);
            }
            indptr.push(indices.len());
        }
        Self::try_new(a.nrows, a.ncols, indptr, indices, values).expect("valid axpby CSR")
    }

    /// Sparse matrix product `self · other` (Gustavson's row-by-row algorithm); the result has
    /// sorted columns. Used for Galerkin coarse operators `Pᵀ A P`.
    pub fn matmul(&self, other: &Csr<T>) -> Csr<T> {
        assert_eq!(self.ncols, other.nrows, "inner dimensions differ");
        let n = other.ncols;
        let mut marker = vec![usize::MAX; n];
        let mut acc = vec![T::ZERO; n];
        let mut indptr = Vec::with_capacity(self.nrows + 1);
        let (mut indices, mut values) = (Vec::new(), Vec::new());
        indptr.push(0);
        let mut cols: Vec<u32> = Vec::new();
        for i in 0..self.nrows {
            cols.clear();
            let (ca, va) = self.row(i);
            for (&k, &a) in ca.iter().zip(va) {
                let (cb, vb) = other.row(k as usize);
                for (&j, &b) in cb.iter().zip(vb) {
                    if marker[j as usize] != i {
                        marker[j as usize] = i;
                        acc[j as usize] = a * b;
                        cols.push(j);
                    } else {
                        acc[j as usize] += a * b;
                    }
                }
            }
            cols.sort_unstable();
            for &j in &cols {
                indices.push(j);
                values.push(acc[j as usize]);
            }
            indptr.push(indices.len());
        }
        Self::try_new(self.nrows, n, indptr, indices, values).expect("valid matmul CSR")
    }

    /// Entry-wise map to another scalar type (e.g. widening a real matrix to complex).
    pub fn map<U: Scalar>(&self, f: impl Fn(T) -> U) -> Csr<U> {
        Csr {
            nrows: self.nrows,
            ncols: self.ncols,
            indptr: self.indptr.clone(),
            indices: self.indices.clone(),
            values: self.values.iter().map(|&v| f(v)).collect(),
        }
    }

    pub fn to_complex(&self) -> Csr<C64> {
        self.map(|v| v.to_c64())
    }

    /// `y[k] = (A x)[start + k]` for the rows covered by `yc`.
    #[inline(always)]
    fn rows_into(&self, start: usize, x: &[T], yc: &mut [T]) {
        let end = start + yc.len();
        assert!(end <= self.nrows && x.len() == self.ncols);
        let ip = &self.indptr[start..=end];
        let (idx, val) = (&self.indices[..], &self.values[..]);
        for (yi, w) in yc.iter_mut().zip(ip.windows(2)) {
            let (lo, hi) = (w[0], w[1]);
            let mut s = T::ZERO;
            for p in lo..hi {
                // SAFETY: `try_new` guarantees `indptr` is monotone with `indptr[nrows] == nnz ==
                // indices.len() == values.len()`, so `p < nnz`; every stored column index is
                // `< ncols == x.len()` (asserted above).
                s += unsafe { *val.get_unchecked(p) * *x.get_unchecked(*idx.get_unchecked(p) as usize) };
            }
            *yi = s;
        }
    }

    /// `y ← A x`.
    pub fn spmv(&self, x: &[T], y: &mut [T]) {
        assert_eq!(x.len(), self.ncols);
        assert_eq!(y.len(), self.nrows);
        let par = par::worth_parallel(self.nnz());
        par::for_chunks_mut(y, CHUNK, par, |c, yc| self.rows_into(c * CHUNK, x, yc));
    }

    /// `y ← y + A x`.
    pub fn spmv_acc(&self, x: &[T], y: &mut [T]) {
        assert_eq!(x.len(), self.ncols);
        assert_eq!(y.len(), self.nrows);
        for (i, yi) in y.iter_mut().enumerate() {
            let (c, v) = self.row(i);
            let mut s = T::ZERO;
            for (&j, &a) in c.iter().zip(v) {
                s += a * x[j as usize];
            }
            *yi += s;
        }
    }

    /// `y ← A x`, returning `xᴴ y` (square matrices) with the reduction fused into the product
    /// so `y` is not re-read. Equals `spmv` followed by [`crate::linalg::vecops::dotc`] bit for bit.
    pub fn spmv_dotc(&self, x: &[T], y: &mut [T]) -> T {
        assert_eq!(self.nrows, self.ncols);
        assert_eq!(x.len(), self.ncols);
        assert_eq!(y.len(), self.nrows);
        let n = self.nrows;
        let ys = par::RawCols::new(y, n, 1);
        let parts = par::chunk_map(n.div_ceil(CHUNK), par::worth_parallel(self.nnz()), |c| {
            let (s, e) = (c * CHUNK, ((c + 1) * CHUNK).min(n));
            let yc = ys.chunk_mut(0, s, e);
            self.rows_into(s, x, yc);
            T::k_dotc(&x[s..e], yc)
        });
        parts.into_iter().fold(T::ZERO, |a, b| a + b)
    }
}

impl<T: Scalar> LinearOperator<T> for Csr<T> {
    fn dim(&self) -> usize {
        assert_eq!(self.nrows, self.ncols, "operator must be square");
        self.nrows
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.spmv(x, y)
    }
    fn apply_dotc(&self, x: &[T], y: &mut [T]) -> T {
        self.spmv_dotc(x, y)
    }
}

// ---------------------------------------------------------------------------------------
// Test problems
// ---------------------------------------------------------------------------------------

/// 1-D Laplacian `tridiag(-1, 2, -1)` of order `n`.
pub fn laplacian_1d(n: usize) -> Csr<f64> {
    let mut t = Vec::with_capacity(3 * n);
    for i in 0..n {
        if i > 0 {
            t.push((i, i - 1, -1.0));
        }
        t.push((i, i, 2.0));
        if i + 1 < n {
            t.push((i, i + 1, -1.0));
        }
    }
    Csr::from_triplets(n, n, &t).unwrap()
}

/// 2-D 5-point Laplacian on an `nx × ny` grid (Dirichlet, natural ordering `i + nx·j`).
pub fn laplacian_2d(nx: usize, ny: usize) -> Csr<f64> {
    let n = nx * ny;
    let mut t = Vec::with_capacity(5 * n);
    for j in 0..ny {
        for i in 0..nx {
            let p = i + nx * j;
            if j > 0 {
                t.push((p, p - nx, -1.0));
            }
            if i > 0 {
                t.push((p, p - 1, -1.0));
            }
            t.push((p, p, 4.0));
            if i + 1 < nx {
                t.push((p, p + 1, -1.0));
            }
            if j + 1 < ny {
                t.push((p, p + nx, -1.0));
            }
        }
    }
    Csr::from_triplets(n, n, &t).unwrap()
}

/// 3-D 7-point Laplacian on an `n³` grid.
pub fn laplacian_3d(n: usize) -> Csr<f64> {
    let idx = |i: usize, j: usize, k: usize| i + n * (j + n * k);
    let mut t = Vec::with_capacity(7 * n * n * n);
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let p = idx(i, j, k);
                if k > 0 {
                    t.push((p, idx(i, j, k - 1), -1.0));
                }
                if j > 0 {
                    t.push((p, idx(i, j - 1, k), -1.0));
                }
                if i > 0 {
                    t.push((p, idx(i - 1, j, k), -1.0));
                }
                t.push((p, p, 6.0));
                if i + 1 < n {
                    t.push((p, idx(i + 1, j, k), -1.0));
                }
                if j + 1 < n {
                    t.push((p, idx(i, j + 1, k), -1.0));
                }
                if k + 1 < n {
                    t.push((p, idx(i, j, k + 1), -1.0));
                }
            }
        }
    }
    Csr::from_triplets(n * n * n, n * n * n, &t).unwrap()
}

/// Lumped-style SPD "mass" matrix `diag(2/3·h)` blended with off-diagonals `h/6` (1-D linear
/// finite-element mass matrix scaled by `h`), for generalised-problem tests.
pub fn mass_1d(n: usize, h: f64) -> Csr<f64> {
    let mut t = Vec::with_capacity(3 * n);
    for i in 0..n {
        if i > 0 {
            t.push((i, i - 1, h / 6.0));
        }
        t.push((i, i, 4.0 * h / 6.0));
        if i + 1 < n {
            t.push((i, i + 1, h / 6.0));
        }
    }
    Csr::from_triplets(n, n, &t).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::rng::random_vector;
    use crate::linalg::vecops;

    #[test]
    fn spmv_matches_dense_and_fused_dot_is_bit_identical() {
        let a = laplacian_2d(70, 61).to_complex();
        let n = a.nrows();
        let (x, mut y1, mut y2) = (random_vector::<C64>(n, 1), vec![C64::ZERO; n], vec![C64::ZERO; n]);
        a.spmv(&x, &mut y1);
        let d = a.spmv_dotc(&x, &mut y2);
        assert_eq!(y1, y2);
        assert_eq!(d, vecops::dotc(&x, &y1));
        let dense = a.to_dense();
        for i in (0..n).step_by(311) {
            let want = (0..n).fold(C64::ZERO, |s, j| s + dense[i * n + j] * x[j]);
            assert!((y1[i] - want).norm() < 1e-12);
        }
    }

    #[test]
    fn symmetry_checks_and_axpby() {
        let a = laplacian_1d(30).to_complex();
        assert!(a.is_symmetric(0.0) && a.is_hermitian(0.0));
        let b = a.map(|v| v * C64::new(0.0, 1.0)); // i·A: complex symmetric, not Hermitian
        assert!(b.is_symmetric(0.0) && !b.is_hermitian(1e-12));
        let c = Csr::axpby(C64::ONE, &a, C64::new(0.0, 1.0), &Csr::identity(30));
        assert_eq!(c.nnz(), a.nnz());
        assert_eq!(c.diagonal()[3], C64::new(2.0, 1.0));
        assert_eq!(a.transpose(), a);
    }

    #[test]
    fn matmul_matches_dense() {
        let a = Csr::axpby(C64::ONE, &laplacian_2d(7, 6).to_complex(), C64::new(0.0, 0.5), &Csr::identity(42));
        let b = laplacian_2d(7, 6).to_complex().transpose();
        let c = a.matmul(&b);
        assert!(c.has_sorted_unique_rows());
        let (da, db, dc) = (a.to_dense(), b.to_dense(), c.to_dense());
        for i in 0..42 {
            for j in 0..42 {
                let want = (0..42).fold(C64::ZERO, |s, k| s + da[i * 42 + k] * db[k * 42 + j]);
                assert!((dc[i * 42 + j] - want).norm() < 1e-12);
            }
        }
    }

    #[test]
    fn triplets_sum_duplicates() {
        let m = Csr::from_triplets(2, 2, &[(0, 1, 1.0), (0, 1, 2.0), (1, 0, 4.0)]).unwrap();
        assert_eq!(m.values(), &[3.0, 4.0]);
        assert!(Csr::<f64>::from_triplets(2, 2, &[(0, 2, 1.0)]).is_err());
    }

    #[test]
    fn triplets_accessors_work() {
        let mut t = Coo::with_capacity(2, 3, 2);
        t.push(0, 2, 1.5);
        t.push(1, 0, -2.0);
        assert_eq!(t.nrows(), 2);
        assert_eq!(t.ncols(), 3);
        assert_eq!(t.rows(), &[0, 1]);
        assert_eq!(t.cols(), &[2, 0]);
        assert_eq!(t.vals(), &[1.5, -2.0]);
        let csr = t.to_csr();
        assert_eq!(csr.get(0, 2), 1.5);
        assert_eq!(csr.get(0, 99), 0.0);
    }

    #[test]
    #[should_panic]
    fn triplets_push_rejects_oob() {
        let mut t = Coo::new(2, 2);
        t.push(2, 0, 1.0);
    }

    #[test]
    #[should_panic]
    fn triplets_push_block_rejects_len_mismatch() {
        let mut t = Coo::new(3, 3);
        t.push_block(&[0, 1], &[0, 1], &[1.0, 2.0, 3.0]);
    }

    #[test]
    #[should_panic(expected = "overflows u32")]
    fn from_raw_rejects_u32_truncation() {
        let big = (1usize << 32) + 1;
        let _ = RCsr::from_raw(1, 2, vec![0, 1], vec![big], vec![1.0]);
    }

    #[test]
    #[should_panic(expected = "shape mismatch")]
    fn sum_all_rejects_shape_mismatch() {
        let a = RCsr::zeros(2, 2);
        let b = RCsr::zeros(2, 3);
        let _ = RCsr::sum_all(vec![a, b]);
    }

    #[test]
    #[should_panic(expected = "at least one part")]
    fn sum_all_rejects_empty() {
        let _ = RCsr::sum_all(Vec::<RCsr>::new());
    }
}
