//! Action of the matrix exponential (port item X1; the role of SLEPc `MFN`): `w = exp(tA) v` for a large
//! sparse or matrix-free `A`, by Krylov projection with adaptive sub-stepping (Sidje's EXPOKIT `expv`
//! strategy), and the dense building block [`expm`] (scaling and squaring with the degree-13 Padé
//! approximant, Higham 2005).
//!
//! `A` is only touched through products, so exponential integrators for linear transient problems
//! (`u' = A u + f`: pulses in optomechanics, wave propagation `A = [[0, I], [−K, 0]]`) need no factorisation
//! and no linear solves. Cost per Krylov step: one product and a Gram–Schmidt pass; the basis has `m + 1`
//! vectors. The step `τ` is chosen so that the a-posteriori error estimate
//! `β h_{m+1,m} τ |e_mᵀ φ₁(τH_m) e₁|` stays below `tol·τ` (error per unit time).

use crate::linalg::dense::Dense;
use crate::linalg::linop::LinearOperator;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::{axpy, dotc, norm2};

/// Solves `A X = B` for dense complex matrices by LU with partial pivoting.
fn solve_dense(a: &Dense<C64>, b: &Dense<C64>) -> Option<Dense<C64>> {
    let n = a.rows();
    let mut m = a.clone();
    let mut x = b.clone();
    let scale = m.frobenius().max(f64::MIN_POSITIVE);
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| m[(i, k)].norm().total_cmp(&m[(j, k)].norm()))?;
        if m[(p, k)].norm() <= 1e-300 * scale {
            return None;
        }
        if p != k {
            for c in 0..n {
                let t = m[(k, c)];
                m[(k, c)] = m[(p, c)];
                m[(p, c)] = t;
            }
            for c in 0..x.cols() {
                let t = x[(k, c)];
                x[(k, c)] = x[(p, c)];
                x[(p, c)] = t;
            }
        }
        for i in k + 1..n {
            let f = m[(i, k)] / m[(k, k)];
            for c in k..n {
                let v = m[(k, c)];
                m[(i, c)] -= f * v;
            }
            for c in 0..x.cols() {
                let v = x[(k, c)];
                x[(i, c)] -= f * v;
            }
        }
    }
    for c in 0..x.cols() {
        for i in (0..n).rev() {
            let mut s = x[(i, c)];
            for j in i + 1..n {
                s -= m[(i, j)] * x[(j, c)];
            }
            x[(i, c)] = s / m[(i, i)];
        }
    }
    Some(x)
}

fn norm1(a: &Dense<C64>) -> f64 {
    (0..a.cols()).map(|j| a.col(j).iter().map(|v| v.norm()).sum::<f64>()).fold(0.0, f64::max)
}

fn lin(parts: &[(f64, &Dense<C64>)], n: usize) -> Dense<C64> {
    Dense::from_fn(n, n, |i, j| parts.iter().map(|(c, m)| m[(i, j)].scale(*c)).sum())
}

/// Dense matrix exponential (scaling and squaring, Padé degree 13).
pub fn expm(a: &Dense<C64>) -> Dense<C64> {
    let n = a.rows();
    assert_eq!(a.cols(), n);
    if n == 0 {
        return Dense::zeros(0, 0);
    }
    const B: [f64; 14] = [
        64764752532480000.0,
        32382376266240000.0,
        7771770303897600.0,
        1187353796428800.0,
        129060195264000.0,
        10559470521600.0,
        670442572800.0,
        33522128640.0,
        1323241920.0,
        40840800.0,
        960960.0,
        16380.0,
        182.0,
        1.0,
    ];
    let theta13 = 5.371920351148152;
    let nrm = norm1(a);
    let s = if nrm > theta13 { (nrm / theta13).log2().ceil().max(0.0) as i32 } else { 0 };
    let scaled = Dense::from_fn(n, n, |i, j| a[(i, j)].scale(0.5f64.powi(s)));
    let id = Dense::<C64>::identity(n);
    let a2 = scaled.matmul(&scaled);
    let a4 = a2.matmul(&a2);
    let a6 = a4.matmul(&a2);
    let u_inner = lin(&[(B[13], &a6), (B[11], &a4), (B[9], &a2)], n);
    let u_poly = lin(&[(1.0, &a6.matmul(&u_inner)), (B[7], &a6), (B[5], &a4), (B[3], &a2), (B[1], &id)], n);
    let u = scaled.matmul(&u_poly);
    let v_inner = lin(&[(B[12], &a6), (B[10], &a4), (B[8], &a2)], n);
    let v = lin(&[(1.0, &a6.matmul(&v_inner)), (B[6], &a6), (B[4], &a4), (B[2], &a2), (B[0], &id)], n);
    let num = lin(&[(1.0, &v), (1.0, &u)], n);
    let den = lin(&[(1.0, &v), (-1.0, &u)], n);
    let mut r = solve_dense(&den, &num).expect("Padé denominator is singular");
    for _ in 0..s {
        r = r.matmul(&r);
    }
    r
}

#[derive(Debug, Clone)]
pub struct ExpvOptions {
    /// Krylov dimension per sub-step.
    pub m: usize,
    /// Error tolerance per unit time (absolute, in the 2-norm of the solution update, relative to `‖v‖`).
    pub tol: f64,
    pub max_steps: usize,
}

impl Default for ExpvOptions {
    fn default() -> Self {
        Self { m: 30, tol: 1e-10, max_steps: 10_000 }
    }
}

#[derive(Debug, Clone)]
pub struct ExpvResult<T> {
    pub w: Vec<T>,
    pub steps: usize,
    pub matvecs: usize,
    /// Sum of the local error estimates (relative to `‖v‖`).
    pub error_estimate: f64,
    pub converged: bool,
}

/// `w = exp(t A) v` (`t` real, may be negative).
pub fn expv<T: Scalar, A: LinearOperator<T> + ?Sized>(a: &A, t: f64, v: &[T], opts: &ExpvOptions) -> ExpvResult<T> {
    let n = a.dim();
    assert_eq!(v.len(), n);
    let beta0 = norm2(v);
    let mut w = v.to_vec();
    if beta0 == 0.0 || t == 0.0 {
        return ExpvResult { w, steps: 0, matvecs: 0, error_estimate: 0.0, converged: true };
    }
    let m = opts.m.min(n.saturating_sub(1)).max(1);
    let sign = t.signum();
    let t_total = t.abs();
    let (mut t_now, mut tau) = (0.0f64, t_total);
    let (mut steps, mut matvecs, mut err_total) = (0usize, 0usize, 0.0f64);
    let mut basis: Vec<Vec<T>> = Vec::with_capacity(m + 1);
    let mut converged = true;
    while t_now < t_total {
        if steps >= opts.max_steps {
            converged = false;
            break;
        }
        // Arnoldi on w
        let beta = norm2(&w);
        if beta == 0.0 {
            break;
        }
        basis.clear();
        basis.push(w.iter().map(|x| x.scale(1.0 / beta)).collect());
        let mut h = Dense::<C64>::zeros(m + 1, m);
        let mut mm = m;
        let mut breakdown = false;
        let mut wv = vec![T::ZERO; n];
        for j in 0..m {
            a.apply(&basis[j], &mut wv);
            matvecs += 1;
            for _ in 0..2 {
                for (i, vi) in basis.iter().enumerate().take(j + 1) {
                    let hij = dotc(vi, &wv);
                    h[(i, j)] += hij.to_c64();
                    axpy(-hij, vi, &mut wv);
                }
            }
            let hn = norm2(&wv);
            h[(j + 1, j)] = C64::new(hn, 0.0);
            if hn <= 1e-14 * beta * (1.0 + h[(0, 0)].norm()) {
                mm = j + 1;
                breakdown = true;
                break;
            }
            basis.push(wv.iter().map(|x| x.scale(1.0 / hn)).collect());
        }
        let hlast = if breakdown { 0.0 } else { h[(m, m - 1)].re };
        // choose τ by the error estimate of φ₁ (augmented exponential)
        let aug = |tau: f64| -> Dense<C64> {
            let k = mm;
            Dense::from_fn(k + 1, k + 1, |i, j| {
                if i < k && j < k {
                    h[(i, j)].scale(sign * tau)
                } else if i == 0 && j == k {
                    C64::ONE
                } else {
                    C64::ZERO
                }
            })
        };
        let mut tau_try = tau.min(t_total - t_now);
        let mut f;
        let mut est;
        let mut shrink = 0;
        loop {
            f = expm(&aug(tau_try));
            // φ₁(τH) e₁ is the last column of the augmented exponential (top mm entries)
            let phi_last = f[(mm - 1, mm)].norm();
            est = if breakdown { 0.0 } else { hlast * tau_try * phi_last };
            if est <= opts.tol * tau_try.max(1e-300) || shrink > 30 {
                break;
            }
            tau_try *= (0.9 * opts.tol * tau_try / est).powf(1.0 / (mm as f64 + 1.0)).clamp(0.2, 0.9);
            shrink += 1;
        }
        // w ← β V exp(τH) e₁
        let mut next = vec![T::ZERO; n];
        for i in 0..mm {
            let c = f[(i, 0)].scale(beta);
            let ct = T::from_parts(c.re, c.im);
            axpy(ct, &basis[i], &mut next);
        }
        w = next;
        err_total += est;
        t_now += tau_try;
        steps += 1;
        // next step: grow if the estimate was comfortable
        tau = if est > 0.0 {
            (tau_try * (opts.tol * tau_try / est).powf(1.0 / (mm as f64 + 1.0)) * 0.9).clamp(tau_try, 10.0 * tau_try)
        } else {
            10.0 * tau_try
        };
        if t_total - t_now < 1e-14 * t_total {
            break;
        }
    }
    ExpvResult { w, steps, matvecs, error_estimate: err_total, converged }
}
