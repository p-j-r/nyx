//! Chunk-level numerical kernels, written once and instantiated for `f64` and `Complex<f64>`.
//!
//! Every kernel works on one *chunk* of rows and uses a fixed accumulator layout, so the result
//! of a reduction depends only on the chunk boundaries (chosen by the caller, see
//! [`crate::linalg::vecops::CHUNK`]) and never on thread count or CPU. This is the same determinism
//! rule as `rsolvers` (`src/kernels.rs`), minus the hand-written SIMD: the portable code below is
//! the reference, and a runtime-dispatched SIMD tier can replace a method for one scalar type
//! later as long as it keeps the summation order (see `docs/STATUS.md`, "Performance").
//!
//! No kernel uses FMA, so results are also identical across CPUs.

use crate::linalg::scalar::{C64, Scalar};

/// Lanes per accumulator group in the dot kernels.
const LANES: usize = 4;

/// The hot loops, overridable per scalar type.
pub trait Kernels: Sized + Copy {
    /// `Σ conj(xᵢ)·yᵢ`.
    fn k_dotc(x: &[Self], y: &[Self]) -> Self;
    /// `Σ xᵢ·yᵢ` (no conjugation; the bilinear form used by COCG).
    fn k_dotu(x: &[Self], y: &[Self]) -> Self;
    /// `y ← y + a·x`.
    fn k_axpy(a: Self, x: &[Self], y: &mut [Self]);
    /// `out[j] = Σᵢ conj(data[j·n + s + i])·w[i]` for `j < out.len()`; `data` holds columns of
    /// length `n` and `w` is the chunk of rows `s..s+w.len()`.
    fn k_dots_cols(data: &[Self], n: usize, s: usize, w: &[Self], out: &mut [Self]);
    /// `w[i] −= Σⱼ data[j·n + s + i]·h[j]`, `j < h.len()`.
    fn k_sub_cols(data: &[Self], n: usize, s: usize, h: &[Self], w: &mut [Self]);
    /// `y[i] += Σⱼ data[j·n + s + i]·c[j]`, `j < c.len()`.
    fn k_add_cols(data: &[Self], n: usize, s: usize, c: &[Self], y: &mut [Self]);
    /// `out[:, c] = src · q[:, c]` for `c < ncols`: a small matrix product on one row tile.
    /// `src` is `len × nrows` and `out` is `len × ncols`, both contiguous column-major; `q` is
    /// `nrows × ncols` with leading dimension `ldq`. This is the restart multiply `V ← V·Q`
    /// (a GEMM), so the real and complex types use `faer`'s matrix product unless the
    /// `portable-kernels` feature is on (see [`crate::linalg::faer_glue`]). The portable form is the
    /// definition; the GEMM may differ in the last bits (FMA, blocking) but, being applied per
    /// fixed row tile on one thread, never depends on the thread count.
    fn k_tile_mult(src: &[Self], len: usize, nrows: usize, q: &[Self], ldq: usize, ncols: usize, out: &mut [Self]);
}

#[cfg_attr(not(any(test, feature = "portable-kernels")), allow(dead_code))]
#[inline(always)]
fn tile_mult_portable<T: Scalar>(
    src: &[T],
    len: usize,
    nrows: usize,
    q: &[T],
    ldq: usize,
    ncols: usize,
    out: &mut [T],
) {
    for c in 0..ncols {
        let o = &mut out[c * len..(c + 1) * len];
        o.fill(T::ZERO);
        combine_generic(src, len, 0, &q[c * ldq..c * ldq + nrows], o, 1.0);
    }
}

#[inline(always)]
#[allow(clippy::chunks_exact_to_as_chunks)] // the lane loop is the measured one (ADR 33); as_chunks would change nothing but the type
fn dot_generic<T: Scalar>(x: &[T], y: &[T], conj_x: bool) -> T {
    debug_assert_eq!(x.len(), y.len());
    let mut acc = [T::ZERO; LANES];
    let (xh, xt) = x.split_at(x.len() / LANES * LANES);
    let (yh, yt) = y.split_at(xh.len());
    if conj_x {
        for (xc, yc) in xh.chunks_exact(LANES).zip(yh.chunks_exact(LANES)) {
            for l in 0..LANES {
                acc[l] += xc[l].conj() * yc[l];
            }
        }
    } else {
        for (xc, yc) in xh.chunks_exact(LANES).zip(yh.chunks_exact(LANES)) {
            for l in 0..LANES {
                acc[l] += xc[l] * yc[l];
            }
        }
    }
    let mut tail = T::ZERO;
    for (a, b) in xt.iter().zip(yt) {
        tail += if conj_x { a.conj() * *b } else { *a * *b };
    }
    (acc[0] + acc[1]) + (acc[2] + acc[3]) + tail
}

#[inline(always)]
fn axpy_generic<T: Scalar>(a: T, x: &[T], y: &mut [T]) {
    debug_assert_eq!(x.len(), y.len());
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += a * *xi;
    }
}

#[inline(always)]
fn dots_cols_generic<T: Scalar>(data: &[T], n: usize, s: usize, w: &[T], out: &mut [T]) {
    for (j, o) in out.iter_mut().enumerate() {
        let col = &data[j * n + s..j * n + s + w.len()];
        *o = dot_generic(col, w, true);
    }
}

/// Combination in column groups of four so `w` is read and written once per group.
#[inline(always)]
fn combine_generic<T: Scalar>(data: &[T], n: usize, s: usize, c: &[T], w: &mut [T], sign: f64) {
    let len = w.len();
    let col = |j: usize| &data[j * n + s..j * n + s + len];
    let mut j = 0;
    while j + 4 <= c.len() {
        let (c0, c1, c2, c3) = (c[j].scale(sign), c[j + 1].scale(sign), c[j + 2].scale(sign), c[j + 3].scale(sign));
        let (v0, v1, v2, v3) = (col(j), col(j + 1), col(j + 2), col(j + 3));
        for i in 0..len {
            // Fixed association order: ((w + c0 v0) + c1 v1) + ... would serialise the dependency
            // chain, so the four products are summed first.
            w[i] += (c0 * v0[i] + c1 * v1[i]) + (c2 * v2[i] + c3 * v3[i]);
        }
        j += 4;
    }
    while j < c.len() {
        axpy_generic(c[j].scale(sign), col(j), w);
        j += 1;
    }
}

macro_rules! impl_kernels {
    ($t:ty) => {
        impl Kernels for $t {
            #[inline]
            fn k_dotc(x: &[Self], y: &[Self]) -> Self {
                dot_generic(x, y, true)
            }
            #[inline]
            fn k_dotu(x: &[Self], y: &[Self]) -> Self {
                dot_generic(x, y, false)
            }
            #[inline]
            fn k_axpy(a: Self, x: &[Self], y: &mut [Self]) {
                axpy_generic(a, x, y)
            }
            #[inline]
            fn k_dots_cols(data: &[Self], n: usize, s: usize, w: &[Self], out: &mut [Self]) {
                dots_cols_generic(data, n, s, w, out)
            }
            #[inline]
            fn k_sub_cols(data: &[Self], n: usize, s: usize, h: &[Self], w: &mut [Self]) {
                combine_generic(data, n, s, h, w, -1.0)
            }
            #[inline]
            fn k_add_cols(data: &[Self], n: usize, s: usize, c: &[Self], y: &mut [Self]) {
                combine_generic(data, n, s, c, y, 1.0)
            }
            #[inline]
            fn k_tile_mult(
                src: &[Self],
                len: usize,
                nrows: usize,
                q: &[Self],
                ldq: usize,
                ncols: usize,
                out: &mut [Self],
            ) {
                #[cfg(not(feature = "portable-kernels"))]
                {
                    crate::linalg::faer_glue::tile_mult(src, len, nrows, q, ldq, ncols, out)
                }
                #[cfg(feature = "portable-kernels")]
                {
                    tile_mult_portable(src, len, nrows, q, ldq, ncols, out)
                }
            }
        }
    };
}

impl_kernels!(f64);
impl_kernels!(C64);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::rng::SplitMix64;

    fn rv<T: Scalar>(n: usize, seed: u64) -> Vec<T> {
        let mut r = SplitMix64::new(seed);
        (0..n).map(|_| T::from_parts(r.next_signed(), r.next_signed())).collect()
    }

    fn check<T: Scalar>() {
        for &n in &[0usize, 1, 3, 4, 5, 17, 1000] {
            let (x, y) = (rv::<T>(n, 1), rv::<T>(n, 2));
            let want_c = x.iter().zip(&y).fold(T::ZERO, |a, (p, q)| a + p.conj() * *q);
            let want_u = x.iter().zip(&y).fold(T::ZERO, |a, (p, q)| a + *p * *q);
            assert!((T::k_dotc(&x, &y) - want_c).modulus() < 1e-12 * (1.0 + want_c.modulus()));
            assert!((T::k_dotu(&x, &y) - want_u).modulus() < 1e-12 * (1.0 + want_u.modulus()));
        }
        // Column kernels against the naive definition, with k = 0..9 to hit the 4-way remainder.
        let (n, s, len) = (50usize, 7usize, 31usize);
        for k in 0..9 {
            let data = rv::<T>(n * 9, 3);
            let w = rv::<T>(len, 4);
            let h = rv::<T>(k, 5);
            let mut out = vec![T::ZERO; k];
            T::k_dots_cols(&data, n, s, &w, &mut out);
            for j in 0..k {
                let want = (0..len).fold(T::ZERO, |a, i| a + data[j * n + s + i].conj() * w[i]);
                assert!((out[j] - want).modulus() < 1e-12);
            }
            let mut sub = w.clone();
            T::k_sub_cols(&data, n, s, &h, &mut sub);
            let mut add = w.clone();
            T::k_add_cols(&data, n, s, &h, &mut add);
            for i in 0..len {
                let comb = (0..k).fold(T::ZERO, |a, j| a + data[j * n + s + i] * h[j]);
                assert!((sub[i] - (w[i] - comb)).modulus() < 1e-12);
                assert!((add[i] - (w[i] + comb)).modulus() < 1e-12);
            }
        }
    }

    fn check_tile_mult<T: Scalar>() {
        let (len, nrows, ncols, ldq) = (37usize, 13usize, 6usize, 15usize);
        let src = rv::<T>(len * nrows, 7);
        let q = rv::<T>(ldq * ncols, 8);
        let (mut got, mut want) = (vec![T::ZERO; len * ncols], vec![T::ZERO; len * ncols]);
        T::k_tile_mult(&src, len, nrows, &q, ldq, ncols, &mut got);
        tile_mult_portable(&src, len, nrows, &q, ldq, ncols, &mut want);
        for (g, w) in got.iter().zip(&want) {
            assert!((*g - *w).modulus() < 1e-12);
        }
        for c in 0..ncols {
            for i in 0..len {
                let naive = (0..nrows).fold(T::ZERO, |a, j| a + src[j * len + i] * q[c * ldq + j]);
                assert!((got[c * len + i] - naive).modulus() < 1e-12);
            }
        }
    }

    #[test]
    fn tile_mult_matches_definition() {
        check_tile_mult::<f64>();
        check_tile_mult::<C64>();
    }

    #[test]
    fn real_kernels() {
        check::<f64>();
    }
    #[test]
    fn complex_kernels() {
        check::<C64>();
    }
}
