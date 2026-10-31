//! Complex Schur decomposition `A = Z T Zᴴ` by Householder Hessenberg reduction and
//! single-shift QR, plus the adjacent-swap reordering used by Krylov–Schur restarts.
//!
//! The matrices involved are small (the Krylov dimension, tens to a few hundred), so clarity and
//! robustness beat speed here: the cost is negligible next to one factorisation or operator
//! application (cf. the 1.7 s of "small eigenproblems" in nyx's 24.85 s elastic example).
//! Results are checked against `faer`'s general eigensolver in the tests.

use super::{Dense, DenseError};
use crate::solvers::scalar::{C64, Scalar};

pub(crate) struct Givens {
    c: f64,
    s: C64,
}

/// Rotation `G = [c s; −s̄ c]` with `G·[a; b] = [r; 0]`.
pub(crate) fn givens(a: C64, b: C64) -> Givens {
    if b == C64::ZERO {
        return Givens { c: 1.0, s: C64::ZERO };
    }
    if a == C64::ZERO {
        return Givens { c: 0.0, s: b.conj().scale(1.0 / b.norm()) };
    }
    let (na, nb) = (a.norm(), b.norm());
    let norm = na.hypot(nb);
    Givens { c: na / norm, s: (a.scale(1.0 / na) * b.conj()).scale(1.0 / norm) }
}

impl Givens {
    /// Rows `r0, r0+1` of `m`, columns `cols`: `[x; y] ← G [x; y]`.
    fn apply_left(&self, m: &mut Dense<C64>, r0: usize, cols: std::ops::Range<usize>) {
        for j in cols {
            let (x, y) = (m[(r0, j)], m[(r0 + 1, j)]);
            m[(r0, j)] = x.scale(self.c) + self.s * y;
            m[(r0 + 1, j)] = y.scale(self.c) - self.s.conj() * x;
        }
    }
    /// Columns `c0, c0+1` of `m`, rows `rows`: `[x y] ← [x y] Gᴴ`.
    fn apply_right(&self, m: &mut Dense<C64>, c0: usize, rows: std::ops::Range<usize>) {
        for i in rows {
            let (x, y) = (m[(i, c0)], m[(i, c0 + 1)]);
            m[(i, c0)] = x.scale(self.c) + self.s.conj() * y;
            m[(i, c0 + 1)] = y.scale(self.c) - self.s * x;
        }
    }
}

/// Reduces `h` to upper Hessenberg form in place, accumulating the transformation in `z`
/// (`z ← z·Q`).
fn hessenberg(h: &mut Dense<C64>, z: &mut Dense<C64>) {
    let n = h.rows();
    for k in 0..n.saturating_sub(2) {
        let x: Vec<C64> = (k + 1..n).map(|i| h[(i, k)]).collect();
        let tail: f64 = x[1..].iter().map(|v| v.norm_sqr()).sum();
        if tail == 0.0 {
            continue;
        }
        let beta = (x[0].norm_sqr() + tail).sqrt();
        let phase = if x[0] == C64::ZERO { C64::ONE } else { x[0].scale(1.0 / x[0].norm()) };
        let mut v = x;
        v[0] += phase.scale(beta); // x − α e₁ with α = −phase·β
        let vn = v.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        for c in &mut v {
            *c = c.scale(1.0 / vn);
        }
        // Left: H ← (I − 2vvᴴ) H on rows k+1.., columns k..
        for j in k..n {
            let w = (0..v.len()).fold(C64::ZERO, |a, i| a + v[i].conj() * h[(k + 1 + i, j)]);
            for i in 0..v.len() {
                h[(k + 1 + i, j)] -= v[i] * w.scale(2.0);
            }
        }
        // Right on H (all rows) and on Z: M ← M (I − 2vvᴴ) on columns k+1..
        for m in [&mut *h, &mut *z] {
            for i in 0..n {
                let w = (0..v.len()).fold(C64::ZERO, |a, j| a + m[(i, k + 1 + j)] * v[j]);
                for j in 0..v.len() {
                    m[(i, k + 1 + j)] -= w.scale(2.0) * v[j].conj();
                }
            }
        }
        for i in k + 2..n {
            h[(i, k)] = C64::ZERO;
        }
    }
}

/// Wilkinson shift from the trailing 2×2 block `[a b; c d]`.
fn wilkinson(a: C64, b: C64, c: C64, d: C64) -> C64 {
    let p = (a - d).scale(0.5);
    let disc = (p * p + b * c).sqrt();
    let denom = if (p.conj() * disc).re >= 0.0 { p + disc } else { p - disc };
    if denom == C64::ZERO { d } else { d - b * c / denom }
}

/// Schur decomposition `A = Z T Zᴴ` with `T` upper triangular and `Z` unitary.
pub fn schur(a: &Dense<C64>) -> Result<(Dense<C64>, Dense<C64>), DenseError> {
    assert_eq!(a.rows(), a.cols());
    if !a.is_finite() {
        return Err(DenseError::NonFinite);
    }
    let n = a.rows();
    let mut h = a.clone();
    let mut z = Dense::identity(n);
    hessenberg(&mut h, &mut z);
    let norm_h = h.frobenius().max(f64::MIN_POSITIVE);
    let eps = f64::EPSILON;
    let mut hi = n.saturating_sub(1);
    let (mut iter, mut total) = (0usize, 0usize);
    while hi > 0 {
        // Find the start `l` of the active unreduced block ending at `hi`.
        let mut l = hi;
        while l > 0 {
            let s = h[(l - 1, l - 1)].norm() + h[(l, l)].norm();
            let s = if s == 0.0 { norm_h } else { s };
            if h[(l, l - 1)].norm() <= eps * s {
                h[(l, l - 1)] = C64::ZERO;
                break;
            }
            l -= 1;
        }
        if l == hi {
            hi -= 1;
            iter = 0;
            continue;
        }
        iter += 1;
        total += 1;
        if iter > 60 || total > 60 * n.max(1) {
            return Err(DenseError::NoConvergence);
        }
        let mu = if iter % 10 == 0 {
            // Exceptional shift to break rare cycles.
            C64::new(h[(hi, hi - 1)].norm() + if hi >= 2 { h[(hi - 1, hi - 2)].norm() } else { 0.0 }, 0.0)
        } else {
            wilkinson(h[(hi - 1, hi - 1)], h[(hi - 1, hi)], h[(hi, hi - 1)], h[(hi, hi)])
        };
        for i in l..=hi {
            h[(i, i)] -= mu;
        }
        let mut rots = Vec::with_capacity(hi - l);
        for k in l..hi {
            let g = givens(h[(k, k)], h[(k + 1, k)]);
            g.apply_left(&mut h, k, k..n);
            h[(k + 1, k)] = C64::ZERO;
            rots.push(g);
        }
        for (off, g) in rots.iter().enumerate() {
            let k = l + off;
            g.apply_right(&mut h, k, 0..k + 2);
            g.apply_right(&mut z, k, 0..n);
        }
        for i in l..=hi {
            h[(i, i)] += mu;
        }
    }
    for j in 0..n {
        for i in j + 1..n {
            h[(i, j)] = C64::ZERO;
        }
    }
    Ok((h, z))
}

/// Swaps the adjacent diagonal entries `t[q,q]` and `t[q+1,q+1]` by a unitary similarity,
/// updating `z` so that `A = Z T Zᴴ` still holds.
fn swap_adjacent(t: &mut Dense<C64>, z: &mut Dense<C64>, q: usize) {
    let n = t.rows();
    let (a, b, c) = (t[(q, q)], t[(q + 1, q + 1)], t[(q, q + 1)]);
    // [c; b−a] is the eigenvector of the 2×2 block for the eigenvalue b.
    let g = givens(c, b - a);
    g.apply_left(t, q, q..n);
    g.apply_right(t, q, 0..q + 2);
    g.apply_right(z, q, 0..n);
    t[(q + 1, q)] = C64::ZERO;
}

/// Reorders the Schur form so the diagonal follows `key` (ascending: smaller key first).
/// Operates on the leading `m × m` block of `t` (with `z` carrying the same `m` columns), and
/// returns the permutation applied (`perm[p]` = original position of the entry now at `p`).
pub fn reorder(t: &mut Dense<C64>, z: &mut Dense<C64>, m: usize, key: impl Fn(C64) -> f64) -> Vec<usize> {
    let mut ids: Vec<usize> = (0..m).collect();
    for p in 0..m {
        let mut best = p;
        for j in p + 1..m {
            if key(t[(j, j)]) < key(t[(best, best)]) {
                best = j;
            }
        }
        for q in (p..best).rev() {
            swap_adjacent(t, z, q);
            ids.swap(q, q + 1);
        }
    }
    ids
}

/// Eigenvector `y` (length `m`, `y[i] = 1`, zero below `i`) of the upper triangular `t` for the
/// eigenvalue `t[i,i]`. Near-repeated eigenvalues are regularised with a tiny perturbation of
/// the denominator, the standard treatment in `trevc`.
pub fn triangular_eigenvector(t: &Dense<C64>, m: usize, i: usize) -> Vec<C64> {
    let mut y = vec![C64::ZERO; m];
    y[i] = C64::ONE;
    let lam = t[(i, i)];
    let smin = f64::EPSILON * t.frobenius().max(f64::MIN_POSITIVE);
    for j in (0..i).rev() {
        let mut rhs = -t[(j, i)];
        for k in j + 1..i {
            rhs -= t[(j, k)] * y[k];
        }
        let mut d = t[(j, j)] - lam;
        if d.norm() < smin {
            d = C64::new(smin, 0.0);
        }
        y[j] = rhs / d;
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solvers::faer_glue::FaerScalar;
    use crate::solvers::rng::SplitMix64;

    fn random(n: usize, seed: u64) -> Dense<C64> {
        let mut r = SplitMix64::new(seed);
        Dense::from_fn(n, n, |_, _| C64::new(r.next_signed(), r.next_signed()))
    }

    fn residual(a: &Dense<C64>, t: &Dense<C64>, z: &Dense<C64>) -> f64 {
        let rec = z.matmul(t).matmul(&z.adjoint());
        let mut d = rec;
        for j in 0..a.cols() {
            for i in 0..a.rows() {
                d[(i, j)] -= a[(i, j)];
            }
        }
        d.frobenius() / a.frobenius().max(1e-300)
    }

    fn unitarity(z: &Dense<C64>) -> f64 {
        let mut p = z.adjoint().matmul(z);
        for i in 0..z.rows() {
            p[(i, i)] -= C64::ONE;
        }
        p.frobenius()
    }

    #[test]
    fn schur_of_random_matrices() {
        for &n in &[1usize, 2, 3, 8, 31, 70] {
            let a = random(n, n as u64);
            let (t, z) = schur(&a).unwrap();
            assert!(residual(&a, &t, &z) < 1e-13, "n={n}");
            assert!(unitarity(&z) < 1e-13);
            for j in 0..n {
                for i in j + 1..n {
                    assert_eq!(t[(i, j)], C64::ZERO);
                }
            }
        }
    }

    #[test]
    fn schur_eigenvalues_match_faer() {
        let n = 25;
        let a = random(n, 7);
        let (t, _) = schur(&a).unwrap();
        let fa = faer::Mat::<faer::c64>::from_fn(n, n, |i, j| a[(i, j)].to_faer());
        let mut want: Vec<C64> = fa.eigenvalues().unwrap();
        let mut got: Vec<C64> = (0..n).map(|i| t[(i, i)]).collect();
        let key = |c: &C64| (c.re, c.im);
        want.sort_by(|x, y| key(x).partial_cmp(&key(y)).unwrap());
        got.sort_by(|x, y| key(x).partial_cmp(&key(y)).unwrap());
        for (g, w) in got.iter().zip(&want) {
            assert!((g - w).norm() < 1e-11, "{g} vs {w}");
        }
    }

    #[test]
    fn schur_special_matrices() {
        // zero, identity, a Jordan block (defective), a diagonal matrix with repeated entries,
        // and a pure permutation (all eigenvalues on the unit circle).
        let zero = Dense::<C64>::zeros(5, 5);
        let jordan = Dense::from_fn(6, 6, |i, j| {
            if i == j {
                C64::new(2.0, 1.0)
            } else if j == i + 1 {
                C64::ONE
            } else {
                C64::ZERO
            }
        });
        let diag = Dense::from_fn(6, 6, |i, j| if i == j { C64::new((i / 2) as f64, 0.0) } else { C64::ZERO });
        let perm = Dense::from_fn(7, 7, |i, j| if j == (i + 1) % 7 { C64::ONE } else { C64::ZERO });
        for a in [Dense::identity(4), zero, jordan, diag, perm] {
            let (t, z) = schur(&a).unwrap();
            assert!(residual(&a, &t, &z) < 1e-12);
            assert!(unitarity(&z) < 1e-13);
        }
    }

    #[test]
    fn reorder_sorts_diagonal_and_keeps_decomposition() {
        let n = 30;
        let a = random(n, 11);
        let (mut t, mut z) = schur(&a).unwrap();
        let mut sorted: Vec<f64> = (0..n).map(|i| -t[(i, i)].re).collect();
        sorted.sort_by(f64::total_cmp);
        let perm = reorder(&mut t, &mut z, n, |c| -c.re);
        for i in 0..n {
            assert!((-t[(i, i)].re - sorted[i]).abs() < 1e-12);
        }
        assert!(residual(&a, &t, &z) < 1e-12);
        assert!(unitarity(&z) < 1e-12);
        let mut seen = perm.clone();
        seen.sort();
        assert_eq!(seen, (0..n).collect::<Vec<_>>());
        // Partial reorder of a leading block keeps the trailing part intact as a decomposition.
        let (mut t2, mut z2) = schur(&a).unwrap();
        reorder(&mut t2, &mut z2, 12, |c| c.im);
        assert!(residual(&a, &t2, &z2) < 1e-12);
    }

    #[test]
    fn triangular_eigenvectors() {
        let n = 20;
        let a = random(n, 3);
        let (t, z) = schur(&a).unwrap();
        for i in 0..n {
            let y = triangular_eigenvector(&t, n, i);
            // x = Z y is an eigenvector of A.
            let x: Vec<C64> = (0..n).map(|r| (0..n).fold(C64::ZERO, |s, k| s + z[(r, k)] * y[k])).collect();
            let ax: Vec<C64> = (0..n).map(|r| (0..n).fold(C64::ZERO, |s, k| s + a[(r, k)] * x[k])).collect();
            let res: f64 = ax.iter().zip(&x).map(|(p, q)| (*p - t[(i, i)] * *q).norm_sqr()).sum::<f64>().sqrt();
            let nx: f64 = x.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
            assert!(res < 1e-10 * nx * a.frobenius(), "i={i} res={res}");
        }
    }

    #[test]
    fn hermitian_eig_matches_definition() {
        let n = 24;
        let r = random(n, 5);
        let h = Dense::from_fn(n, n, |i, j| (r[(i, j)] + r[(j, i)].conj()).scale(0.5));
        let (w, q) = crate::solvers::dense::hermitian_eig(&h).unwrap();
        assert!(w.windows(2).all(|p| p[0] <= p[1]));
        let mut d = Dense::zeros(n, n);
        for i in 0..n {
            d[(i, i)] = C64::new(w[i], 0.0);
        }
        let rec = q.matmul(&d).matmul(&q.adjoint());
        let mut diff = rec;
        for j in 0..n {
            for i in 0..n {
                diff[(i, j)] -= h[(i, j)];
            }
        }
        assert!(diff.frobenius() < 1e-12 * h.frobenius());
        // The same for a real symmetric matrix.
        let rs = Dense::from_fn(n, n, |i, j| ((i + j) as f64).sin() + if i == j { 3.0 } else { 0.0 });
        let (w, q) = crate::solvers::dense::hermitian_eig(&rs).unwrap();
        for (k, wk) in w.iter().enumerate() {
            let col = q.col(k);
            for i in (0..n).step_by(5) {
                let aq: f64 = (0..n).map(|j| rs[(i, j)] * col[j]).sum();
                assert!((aq - wk * col[i]).abs() < 1e-11);
            }
        }
    }
}
