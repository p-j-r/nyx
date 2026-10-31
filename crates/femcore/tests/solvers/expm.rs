//! Matrix exponential and its action: dense identities, modal closed forms for heat, Schrödinger and
//! wave equations, matrix-free operators and the time direction.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::solvers::csr::laplacian_1d;
use femcore::solvers::dense::Dense;
use femcore::solvers::expm::{ExpvOptions, expm, expv};
use femcore::solvers::linop::FnOperator;
use femcore::solvers::{C64, Csr};
use std::f64::consts::PI;

fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

fn sine(n: usize, k: usize) -> Vec<f64> {
    (1..=n).map(|j| (k as f64 * PI * j as f64 / (n as f64 + 1.0)).sin()).collect()
}

fn max_err(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

#[test]
fn dense_expm_closed_forms_and_inverse_property() {
    // nilpotent
    let a = Dense::from_fn(2, 2, |i, j| if (i, j) == (0, 1) { C64::ONE } else { C64::ZERO });
    let e = expm(&a);
    assert!((e[(0, 1)] - C64::ONE).norm() < 1e-14 && (e[(0, 0)] - C64::ONE).norm() < 1e-14 && e[(1, 0)].norm() < 1e-14);
    // rotation generator with a large norm (forces scaling and squaring): exp([[0,−θ],[θ,0]])
    let th = 37.0;
    let r = Dense::from_fn(2, 2, |i, j| match (i, j) {
        (0, 1) => c(-th, 0.0),
        (1, 0) => c(th, 0.0),
        _ => C64::ZERO,
    });
    let e = expm(&r);
    assert!(
        (e[(0, 0)] - c(th.cos(), 0.0)).norm() < 1e-12 && (e[(1, 0)] - c(th.sin(), 0.0)).norm() < 1e-12,
        "{:?}",
        (e[(0, 0)], e[(1, 0)])
    );
    // exp(A) exp(−A) = I for a non-normal complex matrix of norm ~ 30
    let a = Dense::from_fn(6, 6, |i, j| c(((i * 7 + j * 3) % 11) as f64 - 5.0, 0.3 * ((i + 2 * j) % 5) as f64));
    let ma = Dense::from_fn(6, 6, |i, j| -a[(i, j)]);
    let p = expm(&a).matmul(&expm(&ma));
    let mut worst = 0.0f64;
    for i in 0..6 {
        for j in 0..6 {
            worst = worst.max((p[(i, j)] - if i == j { C64::ONE } else { C64::ZERO }).norm());
        }
    }
    assert!(worst < 1e-8, "{worst}");
    // diagonal
    let d = Dense::from_fn(3, 3, |i, j| if i == j { c(-1.0 + i as f64, 0.5 * i as f64) } else { C64::ZERO });
    let e = expm(&d);
    for i in 0..3 {
        let z = d[(i, i)];
        let want = C64::new(z.re.exp() * z.im.cos(), z.re.exp() * z.im.sin());
        assert!((e[(i, i)] - want).norm() < 1e-13);
    }
}

#[test]
fn heat_equation_action_matches_the_modal_solution() {
    // a random start vector excites all 400 modes, so the Krylov space is genuinely large
    let n = 400;
    let l = laplacian_1d(n).map(|v| -v); // A = −L
    let ev = laplacian_1d_eigenvalues(n);
    let v = femcore::solvers::rng::random_vector::<f64>(n, 17);
    let modes: Vec<Vec<f64>> = (1..=n).map(|k| sine(n, k)).collect();
    let coef: Vec<f64> = modes
        .iter()
        .map(|s| s.iter().zip(&v).map(|(a, b)| a * b).sum::<f64>() / s.iter().map(|a| a * a).sum::<f64>())
        .collect();
    for t in [0.01, 0.5, 5.0, 300.0] {
        let r = expv(&l, t, &v, &ExpvOptions::default());
        assert!(r.converged);
        let mut want = vec![0.0; n];
        for k in 0..n {
            let w = coef[k] * (-ev[k] * t).exp();
            for i in 0..n {
                want[i] += w * modes[k][i];
            }
        }
        let e = max_err(&r.w, &want);
        eprintln!("heat t = {t}: error {e:.2e}, {} sub-steps, {} matvecs", r.steps, r.matvecs);
        assert!(e < 1e-8, "t = {t}: {e}");
    }
}

#[test]
fn schrodinger_evolution_is_unitary_and_exact() {
    let n = 300;
    let h = laplacian_1d(n).to_complex();
    let ia = h.map(|z| c(0.0, -1.0) * z); // A = −iH, exp(tA) unitary
    let ev = laplacian_1d_eigenvalues(n);
    let (ka, kb) = (5usize, 60usize);
    let v: Vec<C64> = (0..n).map(|i| c(sine(n, ka)[i], 0.2 * sine(n, kb)[i])).collect();
    let t = 40.0;
    let r = expv(&ia, t, &v, &ExpvOptions { tol: 1e-11, ..Default::default() });
    assert!(r.converged);
    assert!((norm(&r.w) / norm(&v) - 1.0).abs() < 1e-8, "norm drift {}", norm(&r.w) / norm(&v) - 1.0);
    let ph = |lam: f64| c((lam * t).cos(), -(lam * t).sin());
    let want: Vec<C64> = (0..n)
        .map(|i| ph(ev[ka - 1]) * c(sine(n, ka)[i], 0.0) + ph(ev[kb - 1]) * c(0.0, 0.2 * sine(n, kb)[i]))
        .collect();
    let e = r.w.iter().zip(&want).map(|(a, b)| (*a - *b).norm()).fold(0.0, f64::max);
    assert!(e < 1e-7, "{e}");
}

#[test]
fn matrix_free_wave_equation_and_negative_time() {
    // u'' = −L u as a first-order system z = (u, u'), A z = (u', −L u): matrix-free operator
    let n = 200;
    let l = laplacian_1d(n);
    let ev = laplacian_1d_eigenvalues(n);
    let k = 7usize;
    let om = ev[k - 1].sqrt();
    let op = FnOperator::new(2 * n, |z: &[f64], y: &mut [f64]| {
        y[..n].copy_from_slice(&z[n..]);
        l.spmv(&z[..n], &mut y[n..]);
        for v in &mut y[n..] {
            *v = -*v;
        }
    });
    let mut z0 = vec![0.0; 2 * n];
    z0[..n].copy_from_slice(&sine(n, k));
    let t = 12.0;
    let fwd = expv(&op, t, &z0, &ExpvOptions::default());
    for i in 0..n {
        let s = sine(n, k)[i];
        assert!((fwd.w[i] - s * (om * t).cos()).abs() < 1e-8);
        assert!((fwd.w[n + i] + om * s * (om * t).sin()).abs() < 1e-8);
    }
    // backward in time returns to the start
    let back = expv(&op, -t, &fwd.w, &ExpvOptions::default());
    assert!(max_err(&back.w, &z0) < 1e-7, "{}", max_err(&back.w, &z0));
    // sign handling: exp(−tL)v == exp(t(−L))v
    let v = sine(n, 4);
    let a1 = l.clone();
    let a2 = l.map(|x| -x);
    let r1 = expv(&a1, -2.0, &v, &ExpvOptions::default());
    let r2 = expv(&a2, 2.0, &v, &ExpvOptions::default());
    assert!(max_err(&r1.w, &r2.w) < 1e-9);
    // zero time and zero vector
    assert_eq!(expv(&a1, 0.0, &v, &ExpvOptions::default()).w, v);
    assert!(expv(&a1, 1.0, &vec![0.0; n], &ExpvOptions::default()).w.iter().all(|x| *x == 0.0));
    let _ = Csr::<f64>::identity(1);
    let _ = C64::ZERO;
}
