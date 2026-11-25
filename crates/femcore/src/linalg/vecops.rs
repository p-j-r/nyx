//! Level-1 operations on vectors of any [`Scalar`], with the fixed-chunk reduction order.

use crate::linalg::par;
use crate::linalg::scalar::Scalar;

/// Rows per reduction chunk. Part of the numerical contract: changing it changes the last bits
/// of every dot product.
pub const CHUNK: usize = 2048;

/// `Σ conj(xᵢ)·yᵢ`.
pub fn dotc<T: Scalar>(x: &[T], y: &[T]) -> T {
    assert_eq!(x.len(), y.len());
    let nch = x.len().div_ceil(CHUNK);
    let parts = par::chunk_map(nch, par::worth_parallel(x.len()), |c| {
        let (s, e) = (c * CHUNK, ((c + 1) * CHUNK).min(x.len()));
        T::k_dotc(&x[s..e], &y[s..e])
    });
    parts.into_iter().fold(T::ZERO, |a, b| a + b)
}

/// `Σ xᵢ·yᵢ` without conjugation.
pub fn dotu<T: Scalar>(x: &[T], y: &[T]) -> T {
    assert_eq!(x.len(), y.len());
    let nch = x.len().div_ceil(CHUNK);
    let parts = par::chunk_map(nch, par::worth_parallel(x.len()), |c| {
        let (s, e) = (c * CHUNK, ((c + 1) * CHUNK).min(x.len()));
        T::k_dotu(&x[s..e], &y[s..e])
    });
    parts.into_iter().fold(T::ZERO, |a, b| a + b)
}

/// Squared 2-norm `Σ |xᵢ|²`.
pub fn norm2_sq<T: Scalar>(x: &[T]) -> f64 {
    dotc(x, x).re()
}

/// 2-norm (no over/underflow protection beyond `f64` range).
pub fn norm2<T: Scalar>(x: &[T]) -> f64 {
    norm2_sq(x).sqrt()
}

/// `y ← y + a·x`.
pub fn axpy<T: Scalar>(a: T, x: &[T], y: &mut [T]) {
    assert_eq!(x.len(), y.len());
    let par = par::worth_parallel(x.len());
    par::for_chunks_mut(y, CHUNK, par, |c, yc| {
        let s = c * CHUNK;
        T::k_axpy(a, &x[s..s + yc.len()], yc);
    });
}

/// `y ← a·y`.
pub fn scale<T: Scalar>(a: T, y: &mut [T]) {
    let par = par::worth_parallel(y.len());
    par::for_chunks_mut(y, CHUNK, par, |_, yc| {
        for v in yc {
            *v *= a;
        }
    });
}

/// `y ← r·y` for a real `r`.
pub fn scale_real<T: Scalar>(r: f64, y: &mut [T]) {
    let par = par::worth_parallel(y.len());
    par::for_chunks_mut(y, CHUNK, par, |_, yc| {
        for v in yc {
            *v = v.scale(r);
        }
    });
}

// ---------------------------------------------------------------------------------------
// Fused kernels for the Krylov linear solvers
// ---------------------------------------------------------------------------------------

/// `y ← x + β·y` (PETSc `VecAYPX`).
pub fn aypx<T: Scalar>(beta: T, x: &[T], y: &mut [T]) {
    assert_eq!(x.len(), y.len());
    let par = par::worth_parallel(x.len());
    par::for_chunks_mut(y, CHUNK, par, |c, yc| {
        let xs = &x[c * CHUNK..c * CHUNK + yc.len()];
        for (yi, xi) in yc.iter_mut().zip(xs) {
            *yi = *xi + beta * *yi;
        }
    });
}

/// `(xᴴy or xᵀy, ‖x‖²)` in one pass.
pub fn dot_and_norm_sq<T: Scalar>(x: &[T], y: &[T], conjugate: bool) -> (T, f64) {
    assert_eq!(x.len(), y.len());
    let parts = par::chunk_map(x.len().div_ceil(CHUNK), par::worth_parallel(x.len()), |c| {
        let (s, e) = (c * CHUNK, ((c + 1) * CHUNK).min(x.len()));
        let d = if conjugate { T::k_dotc(&x[s..e], &y[s..e]) } else { T::k_dotu(&x[s..e], &y[s..e]) };
        (d, T::k_dotc(&x[s..e], &x[s..e]).re())
    });
    parts.into_iter().fold((T::ZERO, 0.0), |a, b| (a.0 + b.0, a.1 + b.1))
}

/// One fused pass of the (preconditioned) CG update, for the identity (`diag = None`, `z`
/// untouched) or a diagonal preconditioner (`diag = Some(d)`, `z = d ⊙ r`):
/// `x += αp; r −= αw; z = d ⊙ r`, returning `(zᴴr or zᵀr, ‖z‖², ‖r‖²)` of the new vectors.
#[allow(clippy::too_many_arguments)]
pub fn cg_update<T: Scalar>(
    alpha: T,
    p: &[T],
    x: &mut [T],
    w: &[T],
    r: &mut [T],
    z: &mut [T],
    diag: Option<&[T]>,
    conjugate: bool,
) -> (T, f64, f64) {
    let n = x.len();
    assert!(
        p.len() == n && w.len() == n && r.len() == n && (diag.is_none() || (z.len() == n && diag.unwrap().len() == n))
    );
    let (xs, rs, zs) =
        (par::RawCols::new(x, n, 1), par::RawCols::new(r, n, 1), par::RawCols::new(z, n, usize::from(diag.is_some())));
    let parts = par::chunk_map(n.div_ceil(CHUNK), par::worth_parallel(n), |c| {
        let (s, e) = (c * CHUNK, ((c + 1) * CHUNK).min(n));
        let (xc, rc) = (xs.chunk_mut(0, s, e), rs.chunk_mut(0, s, e));
        T::k_axpy(alpha, &p[s..e], xc);
        T::k_axpy(-alpha, &w[s..e], rc);
        match diag {
            Some(d) => {
                let zc = zs.chunk_mut(0, s, e);
                for ((zi, di), ri) in zc.iter_mut().zip(&d[s..e]).zip(rc.iter()) {
                    *zi = *di * *ri;
                }
                let zr = if conjugate { T::k_dotc(zc, rc) } else { T::k_dotu(zc, rc) };
                (zr, T::k_dotc(zc, zc).re(), T::k_dotc(rc, rc).re())
            }
            None => {
                let rr = T::k_dotc(rc, rc).re();
                (if conjugate { T::from_f64(rr) } else { T::k_dotu(rc, rc) }, rr, rr)
            }
        }
    });
    parts.into_iter().fold((T::ZERO, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::rng::random_vector;
    use crate::linalg::scalar::C64;

    fn check<T: Scalar>() {
        let n = 3 * CHUNK + 17;
        let (x, y) = (random_vector::<T>(n, 1), random_vector::<T>(n, 2));
        let want = x.iter().zip(&y).fold(T::ZERO, |a, (p, q)| a + p.conj() * *q);
        assert!((dotc(&x, &y) - want).modulus() < 1e-10);
        let nrm = norm2(&x);
        assert!((nrm * nrm - x.iter().map(|v| v.modulus_sq()).sum::<f64>()).abs() < 1e-8);
        let mut z = y.clone();
        axpy(T::from_parts(0.5, -0.25), &x, &mut z);
        for i in (0..n).step_by(509) {
            assert!((z[i] - (y[i] + T::from_parts(0.5, -0.25) * x[i])).modulus() < 1e-14);
        }
    }

    #[test]
    fn real_and_complex() {
        check::<f64>();
        check::<C64>();
    }
}
