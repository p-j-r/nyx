//! Alternative sparse storage (port items M1 and M1b):
//!
//! * [`Bsr`]: **block CSR** with square `b × b` blocks (elasticity: `b = 3`; piezoelectric or
//!   `(u, φ)` couplings: `b = 2`, `4`). One column index per block instead of per entry (an `8×` /
//!   `b²×` smaller index array), unrolled fixed-size block products for `b ∈ {1, 2, 3, 4, 6}`, and a
//!   generic fallback. Blocks are stored row-major and may contain explicit zeros.
//! * [`SymCsr`]: **symmetric storage** of the upper triangle, for complex *symmetric* (`Aᵀ = A`, PML) or
//!   Hermitian (`Aᴴ = A`) matrices: about half the memory; the product scatters each stored
//!   off-diagonal entry into two rows, so it runs sequentially (a parallel version needs
//!   colouring or per-thread accumulators and would not be bit-reproducible).
//!
//! Both implement [`LinearOperator`] and convert from/to [`Csr`]. Measured speed (see
//! `examples/spmv_formats.rs` and `docs/LEARNINGS.md`): the gain depends on the memory traffic; the
//! numbers are in the docs and are not the 1.5–2× of the port list on every machine.

use crate::linalg::csr::{Csr, CsrError};
use crate::linalg::linop::LinearOperator;
use crate::linalg::scalar::Scalar;

#[derive(Debug, Clone, PartialEq)]
pub struct Bsr<T: Scalar> {
    nbrows: usize,
    nbcols: usize,
    b: usize,
    indptr: Vec<usize>,
    indices: Vec<u32>,
    /// `indices.len() · b²` values, block after block, each block row-major.
    values: Vec<T>,
}

impl<T: Scalar> Bsr<T> {
    /// Builds block CSR from CSR; dimensions must be multiples of `b`. Every block with at least one
    /// stored entry is kept (zero-filled where the CSR has no entry).
    pub fn from_csr(a: &Csr<T>, b: usize) -> Result<Self, CsrError> {
        if b == 0 || !a.nrows().is_multiple_of(b) || !a.ncols().is_multiple_of(b) {
            return Err(CsrError::BadIndptr("dimensions must be multiples of the block size"));
        }
        let (nbrows, nbcols) = (a.nrows() / b, a.ncols() / b);
        let mut indptr = Vec::with_capacity(nbrows + 1);
        indptr.push(0);
        let (mut indices, mut values): (Vec<u32>, Vec<T>) = (Vec::new(), Vec::new());
        let mut block_of_col = vec![usize::MAX; nbcols];
        for ib in 0..nbrows {
            let first = indices.len();
            // block columns present in this block row, sorted
            let mut cols: Vec<u32> = Vec::new();
            for r in ib * b..(ib + 1) * b {
                for &c in a.row(r).0 {
                    cols.push(c / b as u32);
                }
            }
            cols.sort_unstable();
            cols.dedup();
            for (k, &jb) in cols.iter().enumerate() {
                block_of_col[jb as usize] = first + k;
                indices.push(jb);
            }
            values.resize(indices.len() * b * b, T::ZERO);
            for (rr, r) in (ib * b..(ib + 1) * b).enumerate() {
                let (cs, vs) = a.row(r);
                for (&c, &v) in cs.iter().zip(vs) {
                    let (jb, cc) = (c as usize / b, c as usize % b);
                    let blk = block_of_col[jb];
                    values[blk * b * b + rr * b + cc] = v;
                }
            }
            for &jb in &cols {
                block_of_col[jb as usize] = usize::MAX;
            }
            indptr.push(indices.len());
        }
        Ok(Self { nbrows, nbcols, b, indptr, indices, values })
    }

    /// Back to CSR, dropping exact zeros.
    pub fn to_csr(&self) -> Csr<T> {
        let b = self.b;
        let mut trip = Vec::new();
        for ib in 0..self.nbrows {
            for k in self.indptr[ib]..self.indptr[ib + 1] {
                let jb = self.indices[k] as usize;
                for rr in 0..b {
                    for cc in 0..b {
                        let v = self.values[k * b * b + rr * b + cc];
                        if v != T::ZERO {
                            trip.push((ib * b + rr, jb * b + cc, v));
                        }
                    }
                }
            }
        }
        Csr::from_triplets(self.nbrows * b, self.nbcols * b, &trip).unwrap()
    }

    pub fn block_size(&self) -> usize {
        self.b
    }
    pub fn nrows(&self) -> usize {
        self.nbrows * self.b
    }
    pub fn ncols(&self) -> usize {
        self.nbcols * self.b
    }
    pub fn nnz_blocks(&self) -> usize {
        self.indices.len()
    }
    /// Stored scalars (including explicit zeros inside blocks).
    pub fn stored_values(&self) -> usize {
        self.values.len()
    }
    /// Bytes of values + indices.
    pub fn bytes(&self) -> usize {
        self.values.len() * std::mem::size_of::<T>()
            + self.indices.len() * 4
            + self.indptr.len() * std::mem::size_of::<usize>()
    }

    #[inline(always)]
    #[allow(clippy::chunks_exact_to_as_chunks)] // B is a const generic; the block loop is kept as measured
    fn block_row<const B: usize>(&self, ib: usize, x: &[T], y: &mut [T]) {
        let mut acc = [T::ZERO; B];
        let (lo, hi) = (self.indptr[ib], self.indptr[ib + 1]);
        for (blk, &jb) in self.values[lo * B * B..hi * B * B].chunks_exact(B * B).zip(&self.indices[lo..hi]) {
            let jb = jb as usize;
            let xs: &[T; B] = x[jb * B..jb * B + B].try_into().unwrap();
            for (a, row) in acc.iter_mut().zip(blk.chunks_exact(B)) {
                let row: &[T; B] = row.try_into().unwrap();
                let mut s = T::ZERO;
                for c in 0..B {
                    s += row[c] * xs[c];
                }
                *a += s;
            }
        }
        y[ib * B..ib * B + B].copy_from_slice(&acc);
    }

    fn block_row_generic(&self, ib: usize, x: &[T], y: &mut [T]) {
        let b = self.b;
        for r in 0..b {
            y[ib * b + r] = T::ZERO;
        }
        for k in self.indptr[ib]..self.indptr[ib + 1] {
            let jb = self.indices[k] as usize;
            let blk = &self.values[k * b * b..(k + 1) * b * b];
            for r in 0..b {
                let mut s = T::ZERO;
                for c in 0..b {
                    s += blk[r * b + c] * x[jb * b + c];
                }
                y[ib * b + r] += s;
            }
        }
    }

    /// `y ← A x`.
    pub fn spmv(&self, x: &[T], y: &mut [T]) {
        assert_eq!(x.len(), self.ncols());
        assert_eq!(y.len(), self.nrows());
        match self.b {
            1 => (0..self.nbrows).for_each(|i| self.block_row::<1>(i, x, y)),
            2 => (0..self.nbrows).for_each(|i| self.block_row::<2>(i, x, y)),
            3 => (0..self.nbrows).for_each(|i| self.block_row::<3>(i, x, y)),
            4 => (0..self.nbrows).for_each(|i| self.block_row::<4>(i, x, y)),
            6 => (0..self.nbrows).for_each(|i| self.block_row::<6>(i, x, y)),
            _ => (0..self.nbrows).for_each(|i| self.block_row_generic(i, x, y)),
        }
    }
}

impl<T: Scalar> LinearOperator<T> for Bsr<T> {
    fn dim(&self) -> usize {
        assert_eq!(self.nbrows, self.nbcols, "operator must be square");
        self.nrows()
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.spmv(x, y)
    }
}

/// Upper triangle (with the diagonal) of a symmetric (`Aᵀ = A`) or Hermitian (`Aᴴ = A`) matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct SymCsr<T: Scalar> {
    n: usize,
    hermitian: bool,
    indptr: Vec<usize>,
    indices: Vec<u32>,
    values: Vec<T>,
}

impl<T: Scalar> SymCsr<T> {
    /// Keeps the entries with `col ≥ row`. `tol` is the relative asymmetry allowed (checked with
    /// [`Csr::is_symmetric`] / [`Csr::is_hermitian`]); `hermitian = false` means plain symmetry
    /// (complex symmetric for complex `T`).
    pub fn from_csr(a: &Csr<T>, hermitian: bool, tol: f64) -> Result<Self, CsrError> {
        if a.nrows() != a.ncols() {
            return Err(CsrError::BadIndptr("matrix must be square"));
        }
        let ok = if hermitian { a.is_hermitian(tol) } else { a.is_symmetric(tol) };
        if !ok {
            return Err(CsrError::BadIndptr("matrix is not symmetric/Hermitian to the requested tolerance"));
        }
        let n = a.nrows();
        let (mut indptr, mut indices, mut values) = (vec![0usize], Vec::new(), Vec::new());
        for i in 0..n {
            let (cs, vs) = a.row(i);
            for (&c, &v) in cs.iter().zip(vs) {
                if c as usize >= i {
                    indices.push(c);
                    values.push(v);
                }
            }
            indptr.push(indices.len());
        }
        Ok(Self { n, hermitian, indptr, indices, values })
    }

    pub fn to_csr(&self) -> Csr<T> {
        let mut trip = Vec::with_capacity(2 * self.values.len());
        for i in 0..self.n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                let j = self.indices[k] as usize;
                let v = self.values[k];
                trip.push((i, j, v));
                if j != i {
                    trip.push((j, i, if self.hermitian { v.conj() } else { v }));
                }
            }
        }
        Csr::from_triplets(self.n, self.n, &trip).unwrap()
    }

    pub fn dim(&self) -> usize {
        self.n
    }
    pub fn stored_values(&self) -> usize {
        self.values.len()
    }
    pub fn bytes(&self) -> usize {
        self.values.len() * std::mem::size_of::<T>()
            + self.indices.len() * 4
            + self.indptr.len() * std::mem::size_of::<usize>()
    }

    /// `y ← A x` (sequential; each stored off-diagonal entry contributes to two rows).
    pub fn spmv(&self, x: &[T], y: &mut [T]) {
        assert_eq!(x.len(), self.n);
        assert_eq!(y.len(), self.n);
        y.fill(T::ZERO);
        for i in 0..self.n {
            let mut acc = T::ZERO;
            let xi = x[i];
            for k in self.indptr[i]..self.indptr[i + 1] {
                let j = self.indices[k] as usize;
                let v = self.values[k];
                acc += v * x[j];
                if j != i {
                    y[j] += if self.hermitian { v.conj() * xi } else { v * xi };
                }
            }
            y[i] += acc;
        }
    }
}

impl<T: Scalar> LinearOperator<T> for SymCsr<T> {
    fn dim(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.spmv(x, y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr::{laplacian_2d, laplacian_3d};
    use crate::linalg::scalar::C64;

    fn rand_vec<T: Scalar>(n: usize) -> Vec<T> {
        crate::linalg::rng::random_vector::<T>(n, 7)
    }

    fn close<T: Scalar>(a: &[T], b: &[T]) {
        let scale: f64 = b.iter().map(|v| v.modulus_sq()).sum::<f64>().sqrt();
        let err: f64 = a.iter().zip(b).map(|(x, y)| (*x - *y).modulus_sq()).sum::<f64>().sqrt();
        assert!(err <= 1e-13 * scale.max(1.0), "relative error {}", err / scale);
    }

    /// `A ⊗ B`-like block matrix: every nonzero of `a` becomes a dense `b×b` block.
    pub(crate) fn blockify(a: &Csr<f64>, b: usize) -> Csr<f64> {
        let mut t = Vec::new();
        for i in 0..a.nrows() {
            let (cs, vs) = a.row(i);
            for (&j, &v) in cs.iter().zip(vs) {
                for r in 0..b {
                    for c in 0..b {
                        // symmetric positive block pattern
                        let w = if r == c { 2.0 + r as f64 } else { 0.25 * (1.0 + (r + c) as f64) };
                        t.push((i * b + r, j as usize * b + c, v * w));
                    }
                }
            }
        }
        Csr::from_triplets(a.nrows() * b, a.ncols() * b, &t).unwrap()
    }

    #[test]
    fn bsr_product_and_round_trip_for_all_block_sizes() {
        for b in [1usize, 2, 3, 4, 5, 6] {
            let a = blockify(&laplacian_2d(7, 6), b);
            let bsr = Bsr::from_csr(&a, b).unwrap();
            assert_eq!(bsr.to_csr(), a);
            let x = rand_vec::<f64>(a.ncols());
            let (mut y1, mut y2) = (vec![0.0; a.nrows()], vec![0.0; a.nrows()]);
            a.spmv(&x, &mut y1);
            bsr.spmv(&x, &mut y2);
            close(&y2, &y1);
        }
        assert!(Bsr::from_csr(&laplacian_2d(5, 5), 4).is_err());
    }

    #[test]
    fn bsr_complex_and_block_pattern_with_missing_entries() {
        // a block row with a structurally missing entry inside a block
        let t = vec![
            (0, 0, C64::new(1.0, 1.0)),
            (1, 1, C64::new(2.0, 0.0)),
            (0, 3, C64::new(0.0, 1.0)),
            (3, 0, C64::new(0.0, 1.0)),
            (3, 3, C64::new(4.0, 0.0)),
            (2, 2, C64::ONE),
        ];
        let a = Csr::from_triplets(4, 4, &t).unwrap();
        let bsr = Bsr::from_csr(&a, 2).unwrap();
        assert_eq!(bsr.nnz_blocks(), 4);
        assert_eq!(bsr.to_csr(), a);
        let x = rand_vec::<C64>(4);
        let (mut y1, mut y2) = (vec![C64::ZERO; 4], vec![C64::ZERO; 4]);
        a.spmv(&x, &mut y1);
        bsr.spmv(&x, &mut y2);
        close(&y2, &y1);
    }

    #[test]
    fn symmetric_storage_real_complex_symmetric_and_hermitian() {
        let a = laplacian_3d(6);
        let s = SymCsr::from_csr(&a, false, 1e-14).unwrap();
        assert!(s.stored_values() * 2 < a.nnz() + a.nrows() + 2);
        assert_eq!(s.to_csr(), a);
        let x = rand_vec::<f64>(a.nrows());
        let (mut y1, mut y2) = (vec![0.0; a.nrows()], vec![0.0; a.nrows()]);
        a.spmv(&x, &mut y1);
        s.spmv(&x, &mut y2);
        close(&y2, &y1);
        // complex symmetric: A = L + i diag(d) (+ a complex off-diagonal pair)
        let n = a.nrows();
        let d: Vec<_> = (0..n).map(|i| (i, i, C64::new(0.0, 0.1 * (i % 5) as f64))).collect();
        let extra = [(0, 5, C64::new(0.3, -0.2)), (5, 0, C64::new(0.3, -0.2))];
        let mut all = d;
        all.extend_from_slice(&extra);
        let ac = Csr::axpby(C64::ONE, &a.to_complex(), C64::ONE, &Csr::from_triplets(n, n, &all).unwrap());
        let sc = SymCsr::from_csr(&ac, false, 1e-14).unwrap();
        assert!(SymCsr::from_csr(&ac, true, 1e-14).is_err()); // complex symmetric is not Hermitian
        let xc = rand_vec::<C64>(n);
        let (mut z1, mut z2) = (vec![C64::ZERO; n], vec![C64::ZERO; n]);
        ac.spmv(&xc, &mut z1);
        sc.spmv(&xc, &mut z2);
        close(&z2, &z1);
        // Hermitian
        let herm = [(0, 5, C64::new(0.3, -0.2)), (5, 0, C64::new(0.3, 0.2))];
        let ah = Csr::axpby(C64::ONE, &a.to_complex(), C64::ONE, &Csr::from_triplets(n, n, &herm).unwrap());
        let sh = SymCsr::from_csr(&ah, true, 1e-14).unwrap();
        assert_eq!(sh.to_csr(), ah);
        ah.spmv(&xc, &mut z1);
        sh.spmv(&xc, &mut z2);
        close(&z2, &z1);
    }
}
