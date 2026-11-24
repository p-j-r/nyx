//! Everything that touches `faer`'s dense API.
//!
//! nyx pins its own faer version; if it differs from the one in this crate's `Cargo.toml`, the
//! two places to adapt are this file and `direct.rs` (sparse factorisations). The rest of the
//! crate only sees [`FaerScalar`].

use crate::linalg::dense::{Dense, DenseError};
use crate::linalg::scalar::{C64, Scalar};
use faer::traits::ComplexField;
use faer::{Mat, MatMut, Par};

/// Scalars that `faer` can handle, with the conversions to its element type (identical
/// layouts: `f64` and `Complex<f64>` are faer's own `f64` and `c64`).
pub trait FaerScalar: Scalar {
    type Elem: ComplexField<Real = f64> + Copy;
    fn to_faer(self) -> Self::Elem;
    fn from_faer(x: Self::Elem) -> Self;
    /// Identity reinterpretations of slices (implemented per type, so no `unsafe` is needed).
    fn faer_slice(s: &[Self]) -> &[Self::Elem];
    fn faer_slice_mut(s: &mut [Self]) -> &mut [Self::Elem];
}

impl FaerScalar for f64 {
    type Elem = f64;
    #[inline]
    fn to_faer(self) -> f64 {
        self
    }
    #[inline]
    fn from_faer(x: f64) -> f64 {
        x
    }
    #[inline]
    fn faer_slice(s: &[f64]) -> &[f64] {
        s
    }
    #[inline]
    fn faer_slice_mut(s: &mut [f64]) -> &mut [f64] {
        s
    }
}

impl FaerScalar for C64 {
    type Elem = faer::c64;
    #[inline]
    fn to_faer(self) -> faer::c64 {
        self
    }
    #[inline]
    fn from_faer(x: faer::c64) -> C64 {
        x
    }
    #[inline]
    fn faer_slice(s: &[C64]) -> &[faer::c64] {
        s
    }
    #[inline]
    fn faer_slice_mut(s: &mut [C64]) -> &mut [faer::c64] {
        s
    }
}

pub(crate) fn hermitian_eig<T: FaerScalar>(a: &Dense<T>) -> Result<(Vec<f64>, Dense<T>), DenseError> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::evd::{ComputeEigenvectors, self_adjoint_evd, self_adjoint_evd_scratch};

    assert_eq!(a.rows(), a.cols());
    let n = a.rows();
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    if n == 0 {
        return Ok((Vec::new(), Dense::zeros(0, 0)));
    }
    // Mirror the lower triangle so the input is exactly Hermitian whatever the caller stored.
    let m = Mat::<T::Elem>::from_fn(n, n, |i, j| if i >= j { a[(i, j)].to_faer() } else { a[(j, i)].conj().to_faer() });
    let par = dense_par_dim(n);
    let mut u = Mat::<T::Elem>::zeros(n, n);
    let mut s = faer::diag::Diag::<T::Elem>::zeros(n);
    let mut mem = MemBuffer::new(self_adjoint_evd_scratch::<T::Elem>(
        n,
        ComputeEigenvectors::Yes,
        par,
        Default::default(),
    ));
    let evd = self_adjoint_evd(
        m.as_ref(),
        s.as_mut(),
        Some(u.as_mut()),
        par,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .map_err(|_| DenseError::NoConvergence);
    clear_upper_vector_state();
    evd?;
    let w: Vec<f64> = (0..n).map(|i| Scalar::re(T::from_faer(s[i]))).collect();
    let q = Dense::from_fn(n, n, |i, j| T::from_faer(u[(i, j)]));
    Ok((w, q))
}

/// `out = src · q` on one row tile with `faer`'s matrix product, sequential (one tile is one
/// thread's unit of work). `src` is `len × nrows` and `out` is `len × ncols`, contiguous;
/// `q` is `nrows × ncols` with leading dimension `ldq`.
#[cfg(not(feature = "portable-kernels"))]
pub(crate) fn tile_mult<T: FaerScalar>(
    src: &[T],
    len: usize,
    nrows: usize,
    q: &[T],
    ldq: usize,
    ncols: usize,
    out: &mut [T],
) {
    use faer::linalg::matmul::matmul;
    use faer::{Accum, MatRef};
    let lhs = MatRef::from_column_major_slice(T::faer_slice(&src[..len * nrows]), len, nrows);
    let rhs =
        MatRef::from_column_major_slice_with_stride(T::faer_slice(&q[..ldq * (ncols - 1) + nrows]), nrows, ncols, ldq);
    let dst = MatMut::from_column_major_slice_mut(T::faer_slice_mut(&mut out[..len * ncols]), len, ncols);
    matmul(dst, Accum::Replace, lhs, rhs, <T::Elem as faer::traits::ComplexField>::one_impl(), Par::Seq);
    clear_upper_vector_state();
}

/// Executes `vzeroupper` where available. On the benchmark VM (Cascade Lake class), one small
/// faer GEMM made *all later* baseline-SSE code about 3.6x slower (see
/// `examples/gemm_side_effect.rs`); the cause is dirty upper vector state left by wide kernels
/// without `vzeroupper`. Called after every faer call in this crate; it costs about one cycle.
#[cfg(target_arch = "x86_64")]
#[inline]
pub(crate) fn clear_upper_vector_state() {
    #[target_feature(enable = "avx")]
    fn zeroupper() {
        std::arch::x86_64::_mm256_zeroupper()
    }
    if std::arch::is_x86_feature_detected!("avx") {
        // SAFETY: AVX is available, checked just above; the instruction has no memory effects.
        unsafe { zeroupper() }
    }
}

#[cfg(not(target_arch = "x86_64"))]
#[inline]
pub(crate) fn clear_upper_vector_state() {}

/// Thin SVD of the `n × cols` complex matrix `data` (column-major, `n ≥ cols`): singular values in
/// decreasing order and the left singular vectors as an `n × cols` column-major vector.
pub(crate) fn thin_svd_c64(data: &[C64], n: usize, cols: usize) -> Result<(Vec<f64>, Vec<C64>), DenseError> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::svd::{ComputeSvdVectors, svd, svd_scratch};

    assert!(data.len() >= n * cols && n >= cols);
    if data.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let m = Mat::<faer::c64>::from_fn(n, cols, |i, j| data[j * n + i]);
    let par = dense_par(n * cols);
    let mut s_diag = faer::diag::Diag::<faer::c64>::zeros(cols);
    let mut u = Mat::<faer::c64>::zeros(n, cols);
    let mut mem = MemBuffer::new(svd_scratch::<faer::c64>(
        n,
        cols,
        ComputeSvdVectors::Thin,
        ComputeSvdVectors::No,
        par,
        Default::default(),
    ));
    let res = svd(
        m.as_ref(),
        s_diag.as_mut(),
        Some(u.as_mut()),
        None,
        par,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .map_err(|_| DenseError::NoConvergence);
    clear_upper_vector_state();
    res?;
    let s: Vec<f64> = (0..cols).map(|i| s_diag[i].re).collect();
    let mut out = vec![C64::ZERO; n * cols];
    for j in 0..cols {
        for i in 0..n {
            out[j * n + i] = u[(i, j)];
        }
    }
    Ok((s, out))
}

/// `(singular values, U, V)` of [`thin_svd_uv_c64`], matrices column-major.
pub(crate) type SvdUv = (Vec<f64>, Vec<C64>, Vec<C64>);

/// Thin SVD like [`thin_svd_c64`] but also returning the right singular vectors `V`
/// (`cols × cols`, column-major), so that `A = U diag(s) Vᴴ`.
pub(crate) fn thin_svd_uv_c64(data: &[C64], n: usize, cols: usize) -> Result<SvdUv, DenseError> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::svd::{ComputeSvdVectors, svd, svd_scratch};

    assert!(data.len() >= n * cols && n >= cols);
    if data.iter().any(|v| !v.is_finite()) {
        return Err(DenseError::NonFinite);
    }
    let m = Mat::<faer::c64>::from_fn(n, cols, |i, j| data[j * n + i]);
    let par = dense_par(n * cols);
    let mut s_diag = faer::diag::Diag::<faer::c64>::zeros(cols);
    let mut u = Mat::<faer::c64>::zeros(n, cols);
    let mut v = Mat::<faer::c64>::zeros(cols, cols);
    let mut mem = MemBuffer::new(svd_scratch::<faer::c64>(
        n,
        cols,
        ComputeSvdVectors::Thin,
        ComputeSvdVectors::Thin,
        par,
        Default::default(),
    ));
    let res = svd(
        m.as_ref(),
        s_diag.as_mut(),
        Some(u.as_mut()),
        Some(v.as_mut()),
        par,
        MemStack::new(&mut mem),
        Default::default(),
    )
    .map_err(|_| DenseError::NoConvergence);
    clear_upper_vector_state();
    res?;
    let s: Vec<f64> = (0..cols).map(|i| s_diag[i].re).collect();
    let mut uo = vec![C64::ZERO; n * cols];
    let mut vo = vec![C64::ZERO; cols * cols];
    for j in 0..cols {
        for i in 0..n {
            uo[j * n + i] = u[(i, j)];
        }
        for i in 0..cols {
            vo[j * cols + i] = v[(i, j)];
        }
    }
    Ok((s, uo, vo))
}

/// Generalised eigendecomposition `A x = λ B x` of small dense complex matrices. Eigenvalues with
/// `β = 0` (infinite) are returned as `NaN`.
pub fn generalized_eig_c64(a: &Dense<C64>, b: &Dense<C64>) -> Result<(Vec<C64>, Dense<C64>), DenseError> {
    let n = a.rows();
    assert!(a.cols() == n && b.rows() == n && b.cols() == n);
    if n == 0 {
        return Ok((Vec::new(), Dense::zeros(0, 0)));
    }
    if !a.is_finite() || !b.is_finite() {
        return Err(DenseError::NonFinite);
    }
    if n == 1 {
        // faer 0.24.4's generalised eigensolver panics on a 1×1 pencil (scratch sizing), so the
        // trivial case is handled here.
        let lam = if b[(0, 0)].norm() == 0.0 { C64::new(f64::NAN, f64::NAN) } else { a[(0, 0)] / b[(0, 0)] };
        return Ok((vec![lam], Dense::identity(1)));
    }
    let fa = Mat::<faer::c64>::from_fn(n, n, |i, j| a[(i, j)]);
    let fb = Mat::<faer::c64>::from_fn(n, n, |i, j| b[(i, j)]);
    let g = fa.generalized_eigen(&fb).map_err(|_| DenseError::NoConvergence);
    clear_upper_vector_state();
    let g = g?;
    let (sa, sb, u) = (g.S_a(), g.S_b(), g.U());
    let lam = (0..n)
        .map(|i| {
            let (x, y) = (sa.column_vector()[i], sb.column_vector()[i]);
            if y.norm() == 0.0 { C64::new(f64::NAN, f64::NAN) } else { x / y }
        })
        .collect();
    Ok((lam, Dense::from_fn(n, n, |i, j| u[(i, j)])))
}

// ---- dense kernels of the low-rank multifrontal solver (`blr.rs`) ---------------------------------------

/// Dimension threshold up to which small dense eigensolves, LU, triangular solves, and GEMMs stay
/// strictly `Par::Seq` so they are bit-identical across `NYX_THREADS` (`ADR 11`).
pub(crate) const SMALL_DENSE_SEQ_N: usize = 512;

/// Parallelism for a dense `n × n` decomposition (`n <= 512` is strictly `Par::Seq`).
#[inline]
pub(crate) fn dense_par_dim(n: usize) -> Par {
    if n <= SMALL_DENSE_SEQ_N
        || crate::linalg::direct::deterministic_factor_env()
        || crate::par::threads() <= 1
    {
        Par::Seq
    } else {
        faer::get_global_parallelism()
    }
}

/// Parallelism for dense BLR/SVD work of size `work` (`<= 512 × 512` is strictly `Par::Seq`).
#[inline]
pub(crate) fn dense_par(work: usize) -> Par {
    if work <= SMALL_DENSE_SEQ_N * SMALL_DENSE_SEQ_N
        || crate::linalg::direct::deterministic_factor_env()
        || crate::par::threads() <= 1
    {
        Par::Seq
    } else {
        faer::get_global_parallelism()
    }
}

/// In-place LU with partial pivoting of the `n × n` column-major matrix `a` (unit lower `L` below the diagonal, `U`
/// on and above it). Returns `perm` with `(P A)[i, :] = A[perm[i], :]`, or `None` for an exactly singular or
/// non-finite pivot.
pub(crate) fn dense_lu<T: FaerScalar>(a: &mut [T], n: usize) -> Option<Vec<usize>> {
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::linalg::lu::partial_pivoting::factor::{lu_in_place, lu_in_place_scratch};
    if n == 0 {
        return Some(Vec::new());
    }
    let par = dense_par(n * n);
    let (mut perm, mut inv) = (vec![0usize; n], vec![0usize; n]);
    let mut mem = MemBuffer::new(lu_in_place_scratch::<usize, T::Elem>(n, n, par, Default::default()));
    let m = MatMut::from_column_major_slice_mut(T::faer_slice_mut(a), n, n);
    let _ = lu_in_place(m, &mut perm, &mut inv, par, MemStack::new(&mut mem), Default::default());
    clear_upper_vector_state();
    if (0..n).any(|i| {
        let d = a[i * n + i];
        d == T::ZERO || !d.is_finite()
    }) {
        return None;
    }
    Some(perm)
}

/// `B ← L⁻¹ B` with the unit lower triangle of the packed `n × n` matrix `lu`; `B` is `n × k` column-major.
pub(crate) fn solve_unit_lower<T: FaerScalar>(lu: &[T], n: usize, b: &mut [T], k: usize) {
    if n == 0 || k == 0 {
        return;
    }
    let l = faer::MatRef::from_column_major_slice(T::faer_slice(lu), n, n);
    let bm = MatMut::from_column_major_slice_mut(T::faer_slice_mut(b), n, k);
    faer::linalg::triangular_solve::solve_unit_lower_triangular_in_place(l, bm, dense_par(n * k));
    clear_upper_vector_state();
}

/// `B ← B U⁻¹` with the upper triangle of the packed `n × n` matrix `lu`; `B` is `r × n` column-major.
pub(crate) fn solve_upper_right<T: FaerScalar>(lu: &[T], n: usize, b: &mut [T], r: usize) {
    if n == 0 || r == 0 {
        return;
    }
    // B U⁻¹ = (U⁻ᵀ Bᵀ)ᵀ : lower-triangular solve with the transposed view of U
    let u = faer::MatRef::from_column_major_slice(T::faer_slice(lu), n, n);
    let bm = MatMut::from_column_major_slice_mut(T::faer_slice_mut(b), r, n);
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(
        u.transpose(),
        bm.transpose_mut(),
        dense_par(n * r),
    );
    clear_upper_vector_state();
}

/// `C -= A B` for column-major `A: m × k`, `B: k × n`, `C: m × n` (all contiguous).
pub(crate) fn gemm_sub<T: FaerScalar>(c: &mut [T], a: &[T], b: &[T], m: usize, n: usize, k: usize) {
    use faer::Accum;
    use faer::MatRef;
    use faer::linalg::matmul::matmul;
    if m == 0 || n == 0 || k == 0 {
        return;
    }
    let am = MatRef::from_column_major_slice(T::faer_slice(a), m, k);
    let bm = MatRef::from_column_major_slice(T::faer_slice(b), k, n);
    let cm = MatMut::from_column_major_slice_mut(T::faer_slice_mut(c), m, n);
    matmul(cm, Accum::Add, am, bm, T::to_faer(-T::ONE), dense_par(m * n));
    clear_upper_vector_state();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generalized_eig_handles_tiny_pencils() {
        for n in 1..=5usize {
            let a = Dense::from_fn(n, n, |i, j| {
                C64::new((i + 2 * j) as f64 + if i == j { 3.0 } else { 0.0 }, (i as f64 - j as f64) * 0.1)
            });
            let b =
                Dense::from_fn(n, n, |i, j| if i == j { C64::new(1.0 + i as f64, 0.2) } else { C64::new(0.1, 0.0) });
            let (lam, x) = generalized_eig_c64(&a, &b).unwrap();
            for k in 0..n {
                // A x = λ B x
                for i in 0..n {
                    let ax = (0..n).fold(C64::ZERO, |s, j| s + a[(i, j)] * x[(j, k)]);
                    let bx = (0..n).fold(C64::ZERO, |s, j| s + b[(i, j)] * x[(j, k)]);
                    assert!((ax - lam[k] * bx).norm() < 1e-10 * (1.0 + ax.norm()), "n={n} k={k}");
                }
            }
        }
    }
}
