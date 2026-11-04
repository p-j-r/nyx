//! Column-major block of `n`-vectors with cache-blocked multi-vector kernels (the role of
//! SLEPc's `BV` and PETSc's `VecMDot`/`VecMAXPY`), generic over the scalar and with *block*
//! (several right-hand sides at once) variants for the block Krylov methods.
//!
//! Row loops are blocked in [`BCHUNK`] rows so the vectors being orthogonalised stay cache
//! resident while every basis column streams through once. Reductions use the fixed chunk order
//! described in [`crate::linalg::par`].

use crate::linalg::par::{self, RawCols};
use crate::linalg::scalar::Scalar;

/// Rows per block in the multi-vector kernels (also the reduction chunk).
pub const BCHUNK: usize = 2048;
/// Rows per scratch tile in [`Basis::mult_in_place`].
const TILE_ROWS: usize = 512;

/// `cap` columns of length `n`, stored contiguously column after column.
#[derive(Debug, Clone)]
pub struct Basis<T: Scalar> {
    n: usize,
    cap: usize,
    data: Vec<T>,
}

impl<T: Scalar> Basis<T> {
    pub fn new(n: usize, cap: usize) -> Self {
        Self { n, cap, data: vec![T::ZERO; n * cap] }
    }

    /// Vector length.
    pub fn len(&self) -> usize {
        self.n
    }
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
    /// Number of columns allocated.
    pub fn capacity(&self) -> usize {
        self.cap
    }

    #[inline]
    pub fn col(&self, j: usize) -> &[T] {
        &self.data[j * self.n..(j + 1) * self.n]
    }
    #[inline]
    pub fn col_mut(&mut self, j: usize) -> &mut [T] {
        &mut self.data[j * self.n..(j + 1) * self.n]
    }
    /// Columns `a..b` as one contiguous column-major slice.
    #[inline]
    pub fn cols(&self, a: usize, b: usize) -> &[T] {
        &self.data[a * self.n..b * self.n]
    }
    #[inline]
    pub fn cols_mut(&mut self, a: usize, b: usize) -> &mut [T] {
        &mut self.data[a * self.n..b * self.n]
    }

    pub fn view(&self) -> BasisView<'_, T> {
        BasisView { n: self.n, data: &self.data }
    }

    /// Columns `0..j` read-only, plus all columns from `j` on as one mutable column-major slice
    /// (for orthogonalising freshly computed vectors against their predecessors in place).
    pub fn split_cols_mut(&mut self, j: usize) -> (BasisView<'_, T>, &mut [T]) {
        let n = self.n;
        let (head, tail) = self.data.split_at_mut(j * n);
        (BasisView { n, data: head }, tail)
    }

    /// Column `a` shared and column `b` mutable (`a != b`).
    pub fn col_pair_mut(&mut self, a: usize, b: usize) -> (&[T], &mut [T]) {
        assert_ne!(a, b);
        let n = self.n;
        if a < b {
            let (lo, hi) = self.data.split_at_mut(b * n);
            (&lo[a * n..(a + 1) * n], &mut hi[..n])
        } else {
            let (lo, hi) = self.data.split_at_mut(a * n);
            (&hi[..n], &mut lo[b * n..(b + 1) * n])
        }
    }

    /// Copies column `from` over column `to`.
    pub fn copy_col(&mut self, from: usize, to: usize) {
        if from != to {
            let n = self.n;
            self.data.copy_within(from * n..(from + 1) * n, to * n);
        }
    }

    /// Moves columns `from..from+count` to `to..to+count` (ranges may overlap).
    pub fn move_cols(&mut self, from: usize, to: usize, count: usize) {
        let n = self.n;
        self.data.copy_within(from * n..(from + count) * n, to * n);
    }

    /// `V[:, first..first+ncols_q] ← V[:, first..first+nrows_q] · Q` for a column-major `Q`
    /// with leading dimension `ldq` (`BVMultInPlace`). Done in cache-sized row tiles so the
    /// basis is traversed once and no second basis-sized buffer is needed.
    pub fn mult_in_place(&mut self, first: usize, nrows_q: usize, ncols_q: usize, q: &[T], ldq: usize) {
        assert!(first + nrows_q.max(ncols_q) <= self.cap);
        assert!(ldq >= nrows_q && q.len() >= ldq * ncols_q.max(1) - (ldq - nrows_q));
        if nrows_q == 0 || ncols_q == 0 {
            return;
        }
        let n = self.n;
        let nblocks = n.div_ceil(TILE_ROWS);
        let cols = RawCols::new(&mut self.data, n, self.cap);
        let work = |blk: usize| {
            let s = blk * TILE_ROWS;
            let e = (s + TILE_ROWS).min(n);
            let len = e - s;
            // Gather the source rows into a contiguous scratch tile, multiply, scatter back.
            // Only raw copies touch the shared storage, and every block owns rows `s..e` of all
            // columns, so concurrent blocks never alias.
            let mut src = vec![T::ZERO; len * nrows_q];
            for j in 0..nrows_q {
                // SAFETY: rows s..e of column `first+j` belong to this block only.
                unsafe {
                    std::ptr::copy_nonoverlapping(cols.at(s, first + j), src[j * len..].as_mut_ptr(), len);
                }
            }
            let mut out = vec![T::ZERO; len * ncols_q];
            T::k_tile_mult(&src, len, nrows_q, q, ldq, ncols_q, &mut out);
            for c in 0..ncols_q {
                // SAFETY: as above.
                unsafe {
                    std::ptr::copy_nonoverlapping(out[c * len..].as_ptr(), cols.at(s, first + c), len);
                }
            }
        };
        par::chunk_map(nblocks, par::worth_parallel(n * nrows_q), work);
    }
}

/// Borrowed leading columns of a [`Basis`] with the read-only multi-vector kernels.
#[derive(Clone, Copy)]
pub struct BasisView<'a, T: Scalar> {
    n: usize,
    data: &'a [T],
}

impl<T: Scalar> BasisView<'_, T> {
    /// Vector length.
    pub fn len(&self) -> usize {
        self.n
    }
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
    #[inline]
    pub fn col(&self, j: usize) -> &[T] {
        &self.data[j * self.n..(j + 1) * self.n]
    }
    /// Columns `a..b` as one contiguous column-major slice.
    #[inline]
    pub fn cols(&self, a: usize, b: usize) -> &[T] {
        &self.data[a * self.n..b * self.n]
    }

    /// `h[j + c·k] = V[:, j]ᴴ w[:, c]` for `j < k`, `c < b` (`VecMDot`, block form), where `w`
    /// holds `b` columns of length `n`. If `nrm2` is given, also `nrm2[c] = ‖w[:, c]‖²`, computed
    /// while the rows are in cache.
    pub fn dots_block(&self, k: usize, w: &[T], b: usize, h: &mut [T], mut nrm2: Option<&mut [f64]>) {
        let n = self.n;
        assert!(k * n <= self.data.len() && w.len() >= n * b && h.len() >= k * b);
        let want_nrm = nrm2.is_some();
        let slot = k * b + if want_nrm { b } else { 0 };
        let nch = n.div_ceil(BCHUNK);
        let parts = par::chunk_map(nch, par::worth_parallel(n * (k + b)), |c| {
            let (s, e) = (c * BCHUNK, ((c + 1) * BCHUNK).min(n));
            let mut part = vec![T::ZERO; slot];
            for cc in 0..b {
                let wc = &w[cc * n + s..cc * n + e];
                T::k_dots_cols(self.data, n, s, wc, &mut part[cc * k..(cc + 1) * k]);
                if want_nrm {
                    part[k * b + cc] = T::k_dotc(wc, wc);
                }
            }
            part
        });
        h[..k * b].fill(T::ZERO);
        if let Some(r) = nrm2.as_deref_mut() {
            r[..b].fill(0.0);
        }
        for part in &parts {
            for (o, p) in h[..k * b].iter_mut().zip(&part[..k * b]) {
                *o += *p;
            }
            if let Some(r) = nrm2.as_deref_mut() {
                for cc in 0..b {
                    r[cc] += part[k * b + cc].re();
                }
            }
        }
    }

    /// Single-vector form of [`dots_block`](Self::dots_block).
    pub fn dots(&self, k: usize, w: &[T], h: &mut [T]) {
        self.dots_block(k, w, 1, h, None)
    }

    /// `w[:, c] ← w[:, c] − V[:, :k]·h[:, c]` for `c < b` (`VecMAXPY`, block form).
    pub fn sub_block(&self, k: usize, h: &[T], w: &mut [T], b: usize) {
        let n = self.n;
        assert!(k * n <= self.data.len() && w.len() >= n * b && h.len() >= k * b);
        let nch = n.div_ceil(BCHUNK);
        let wc = RawCols::new(w, n, b);
        par::chunk_map(nch, par::worth_parallel(n * (k + b)), |c| {
            let (s, e) = (c * BCHUNK, ((c + 1) * BCHUNK).min(n));
            for cc in 0..b {
                // SAFETY: chunk `c` owns rows s..e of every column of `w`.
                let dst = unsafe { wc.rows(cc, s, e) };
                T::k_sub_cols(self.data, n, s, &h[cc * k..(cc + 1) * k], dst);
            }
        });
    }

    /// Single-vector form of [`sub_block`](Self::sub_block).
    pub fn sub_combination(&self, k: usize, h: &[T], w: &mut [T]) {
        self.sub_block(k, h, w, 1)
    }

    /// Fused second Gram–Schmidt pass: `w ← w − V h`, then (per row block, while the block of
    /// `V` is still in cache) the product `h2 = Vᴴ w_new` and the squared column norms
    /// `nrm2[c] = ‖w_new[:, c]‖²`.
    ///
    /// Lanczos/Arnoldi vectors lose most of their norm in the first pass, so the refinement
    /// pass is needed on almost every step (≈99% in the `rsolvers` Lanczos measurements);
    /// fusing it saves one sweep over `V`.
    pub fn sub_block_then_dots(&self, k: usize, h: &[T], w: &mut [T], b: usize, h2: &mut [T], nrm2: &mut [f64]) {
        let n = self.n;
        assert!(
            k * n <= self.data.len() && w.len() >= n * b && h.len() >= k * b && h2.len() >= k * b && nrm2.len() >= b
        );
        let nch = n.div_ceil(BCHUNK);
        let slot = k * b + b;
        let wc = RawCols::new(w, n, b);
        let parts = par::chunk_map(nch, par::worth_parallel(n * (k + b)), |c| {
            let (s, e) = (c * BCHUNK, ((c + 1) * BCHUNK).min(n));
            let mut part = vec![T::ZERO; slot];
            for cc in 0..b {
                // SAFETY: chunk `c` owns rows s..e of every column of `w`.
                let dst = unsafe { wc.rows(cc, s, e) };
                T::k_sub_cols(self.data, n, s, &h[cc * k..(cc + 1) * k], dst);
                T::k_dots_cols(self.data, n, s, dst, &mut part[cc * k..(cc + 1) * k]);
                part[k * b + cc] = T::k_dotc(dst, dst);
            }
            part
        });
        h2[..k * b].fill(T::ZERO);
        nrm2[..b].fill(0.0);
        for part in &parts {
            for (o, p) in h2[..k * b].iter_mut().zip(&part[..k * b]) {
                *o += *p;
            }
            for cc in 0..b {
                nrm2[cc] += part[k * b + cc].re();
            }
        }
    }

    /// `y ← y + V[:, :k]·c`.
    pub fn add_combination(&self, k: usize, c: &[T], y: &mut [T]) {
        let n = self.n;
        assert!(k * n <= self.data.len() && y.len() == n && c.len() >= k);
        let nch = n.div_ceil(BCHUNK);
        let yc = RawCols::new(y, n, 1);
        par::chunk_map(nch, par::worth_parallel(n * k), |ch| {
            let (s, e) = (ch * BCHUNK, ((ch + 1) * BCHUNK).min(n));
            // SAFETY: chunk `ch` owns rows s..e of `y`.
            T::k_add_cols(self.data, n, s, &c[..k], unsafe { yc.rows(0, s, e) });
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::rng::random_vector;
    use crate::linalg::scalar::C64;

    fn filled<T: Scalar>(n: usize, k: usize) -> Basis<T> {
        let mut b = Basis::new(n, k);
        for j in 0..k {
            b.col_mut(j).copy_from_slice(&random_vector(n, j as u64 + 1));
        }
        b
    }

    fn check<T: Scalar>() {
        let (n, k, nb) = (5 * BCHUNK + 123, 7, 3);
        let v = filled::<T>(n, k);
        let w = random_vector::<T>(n * nb, 99);
        let mut h = vec![T::ZERO; k * nb];
        let mut nrm = vec![0.0; nb];
        v.view().dots_block(k, &w, nb, &mut h, Some(&mut nrm));
        for c in 0..nb {
            let wc = &w[c * n..(c + 1) * n];
            let want_n: f64 = wc.iter().map(|x| x.modulus_sq()).sum();
            assert!((nrm[c] - want_n).abs() < 1e-8 * want_n);
            for j in 0..k {
                let r = v.col(j).iter().zip(wc).fold(T::ZERO, |a, (p, q)| a + p.conj() * *q);
                assert!((r - h[c * k + j]).modulus() < 1e-9 * (1.0 + r.modulus()));
            }
        }
        // sub_block and its fused variant agree with the naive update, and each other exactly.
        let mut w1 = w.clone();
        v.view().sub_block(k, &h, &mut w1, nb);
        for c in 0..nb {
            for i in (0..n).step_by(997) {
                let comb = (0..k).fold(T::ZERO, |a, j| a + v.col(j)[i] * h[c * k + j]);
                assert!((w1[c * n + i] - (w[c * n + i] - comb)).modulus() < 1e-12);
            }
        }
        let mut w2 = w.clone();
        let (mut h2, mut n2) = (vec![T::ZERO; k * nb], vec![0.0; nb]);
        v.view().sub_block_then_dots(k, &h, &mut w2, nb, &mut h2, &mut n2);
        assert_eq!(w1, w2);
        let (mut h2r, mut n2r) = (vec![T::ZERO; k * nb], vec![0.0; nb]);
        v.view().dots_block(k, &w2, nb, &mut h2r, Some(&mut n2r));
        assert_eq!(h2, h2r);
        assert_eq!(n2, n2r);
        // add_combination undoes sub_combination for one column.
        let mut y = w1[..n].to_vec();
        v.view().add_combination(k, &h[..k], &mut y);
        for i in (0..n).step_by(991) {
            assert!((y[i] - w[i]).modulus() < 1e-12);
        }
    }

    #[test]
    fn kernels_match_naive_real() {
        check::<f64>();
    }
    #[test]
    fn kernels_match_naive_complex() {
        check::<C64>();
    }

    fn mult_check<T: Scalar>() {
        let (n, k) = (3000, 6);
        let mut b = filled::<T>(n, 8);
        let orig = b.clone();
        // Q is 6 x 4 (ldq = 7), applied to columns 1..7 giving columns 1..5.
        let ldq = 7;
        let q = random_vector::<T>(ldq * 4, 5);
        b.mult_in_place(1, k, 4, &q, ldq);
        for jc in 0..4 {
            for i in (0..n).step_by(131) {
                let r = (0..k).fold(T::ZERO, |a, jr| a + orig.col(1 + jr)[i] * q[jc * ldq + jr]);
                assert!((r - b.col(1 + jc)[i]).modulus() < 1e-12);
            }
        }
        assert_eq!(b.col(0), orig.col(0));
        assert_eq!(b.col(7), orig.col(7));
        assert_eq!(b.col(5), orig.col(5));
    }

    #[test]
    fn mult_in_place_matches_naive() {
        mult_check::<f64>();
        mult_check::<C64>();
    }
}
