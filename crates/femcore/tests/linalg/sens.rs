//! Eigenvalue derivatives against finite differences and closed forms, left eigenvectors, and
//! mode tracking through exact and avoided crossings.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, mass_1d};
use femcore::linalg::dense::Dense;
use femcore::linalg::faer_glue::generalized_eig_c64;
use femcore::linalg::nep::{NepTerm, NonlinearEigenProblem};
use femcore::linalg::sens::{Tracker, eigenvalue_derivative, left_eigenvector, match_modes, nonlinear_eigenvalue_derivative};
use femcore::linalg::{C64, Csr, Scalar};

fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

fn dense_of(a: &Csr<C64>) -> Dense<C64> {
    let n = a.nrows();
    let d = a.to_dense();
    Dense::from_fn(n, n, |i, j| d[i * n + j])
}

/// All eigenpairs of `A x = λ B x` by the dense solver (small `n`).
fn eigs(a: &Csr<C64>, b: Option<&Csr<C64>>) -> (Vec<C64>, Vec<Vec<C64>>) {
    let n = a.nrows();
    let bd = b.map_or_else(|| Dense::identity(n), dense_of);
    let (vals, v) = generalized_eig_c64(&dense_of(a), &bd).unwrap();
    (vals, (0..n).map(|j| v.col(j).to_vec()).collect())
}

fn pick(vals: &[C64], target: C64) -> usize {
    (0..vals.len()).min_by(|&i, &j| (vals[i] - target).norm().total_cmp(&(vals[j] - target).norm())).unwrap()
}

#[test]
fn hermitian_generalised_derivative_matches_finite_differences() {
    // A(p) = L + p D (D = diagonal potential), B = FE mass; y = x for the Hermitian pencil
    let n = 40;
    let l = laplacian_1d(n).to_complex();
    let dtrip: Vec<_> = (0..n).map(|i| (i, i, c(1.0 + 0.5 * (i as f64 / n as f64), 0.0))).collect();
    let d = Csr::from_triplets(n, n, &dtrip).unwrap();
    let b = mass_1d(n, 1.0 / (n as f64 + 1.0)).to_complex();
    let a_of = |p: f64| Csr::axpby(C64::ONE, &l, c(p, 0.0), &d);
    let p0 = 0.7;
    let (vals, vecs) = eigs(&a_of(p0), Some(&b));
    for target in [0.5, 3.0, 20.0] {
        let k = pick(&vals, c(target * 40.0, 0.0));
        let (lam, x) = (vals[k], &vecs[k]);
        let der = eigenvalue_derivative(lam, x, x, Some(&b), &d, None);
        let h = 1e-5;
        let (vp, _) = eigs(&a_of(p0 + h), Some(&b));
        let (vm, _) = eigs(&a_of(p0 - h), Some(&b));
        let fd = (vp[pick(&vp, lam)] - vm[pick(&vm, lam)]) / (2.0 * h);
        assert!((der - fd).norm() < 1e-6 * fd.norm().max(1.0), "λ = {lam}: derivative {der} vs finite difference {fd}");
    }
}

#[test]
fn complex_symmetric_derivative_uses_the_conjugate_as_left_vector() {
    // A(p) = L + i p D complex symmetric; y = conj(x) turns yᴴ into xᵀ
    let n = 40;
    let l = laplacian_1d(n).to_complex();
    let dtrip: Vec<_> = (0..n).map(|i| (i, i, c(0.0, 1.0 + (i % 3) as f64))).collect();
    let d = Csr::from_triplets(n, n, &dtrip).unwrap();
    let a_of = |p: f64| Csr::axpby(C64::ONE, &l, c(p, 0.0), &d);
    let p0 = 0.3;
    let (vals, vecs) = eigs(&a_of(p0), None);
    let k = pick(&vals, c(1.0, 1.0));
    let (lam, x) = (vals[k], &vecs[k]);
    let y: Vec<C64> = x.iter().map(|v| v.conj()).collect();
    let der = eigenvalue_derivative(lam, x, &y, None, &d, None);
    let h = 1e-6;
    let (vp, _) = eigs(&a_of(p0 + h), None);
    let (vm, _) = eigs(&a_of(p0 - h), None);
    let fd = (vp[pick(&vp, lam)] - vm[pick(&vm, lam)]) / (2.0 * h);
    assert!((der - fd).norm() < 1e-6 * fd.norm(), "{der} vs {fd}");
}

#[test]
fn non_normal_matrix_needs_a_true_left_eigenvector() {
    // convection–diffusion-like tridiagonal, non-normal: y computed by inverse iteration on Aᴴ
    // (small n: the eigenvalue condition number grows like 5^(n/2), which also limits the finite-difference reference)
    let n = 12;
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, c(2.0 + 0.1 * i as f64, 0.0)));
        if i > 0 {
            t.push((i, i - 1, c(-1.5, 0.0)));
        }
        if i + 1 < n {
            t.push((i, i + 1, c(-0.3, 0.0)));
        }
    }
    let a0 = Csr::from_triplets(n, n, &t).unwrap();
    let dd: Vec<_> = (0..n).map(|i| (i, i, c(((i * 7) % 5) as f64 * 0.1, 0.0))).collect();
    let d = Csr::from_triplets(n, n, &dd).unwrap();
    let a_of = |p: f64| Csr::axpby(C64::ONE, &a0, c(p, 0.0), &d);
    let (vals, vecs) = eigs(&a_of(0.0), None);
    let k = pick(&vals, c(2.5, 0.0));
    let (lam, x) = (vals[k], &vecs[k]);
    let y = left_eigenvector(&a_of(0.0), None, lam, x, 4).unwrap();
    // yᴴ A ≈ λ yᴴ  ⇔  Aᴴ y ≈ λ̄ y
    let ah = a_of(0.0).transpose().map(|v| v.conj());
    let mut ay = vec![C64::ZERO; n];
    ah.spmv(&y, &mut ay);
    let res: f64 = ay.iter().zip(&y).map(|(p, q)| (*p - lam.conj() * *q).norm_sqr()).sum::<f64>().sqrt();
    assert!(res < 1e-8, "left residual {res}");
    let der = eigenvalue_derivative(lam, x, &y, None, &d, None);
    let wrong = eigenvalue_derivative(lam, x, x, None, &d, None); // using the right vector as left one
    let h = 1e-5;
    let (vp, _) = eigs(&a_of(h), None);
    let (vm, _) = eigs(&a_of(-h), None);
    let fd = (vp[pick(&vp, lam)] - vm[pick(&vm, lam)]) / (2.0 * h);
    assert!((der - fd).norm() < 1e-5 * fd.norm().max(1e-3), "{der} vs {fd}");
    assert!(
        (wrong - fd).norm() > 1e-3 * fd.norm().max(1e-3),
        "the right vector should not be a valid left vector here"
    );
}

#[test]
fn quadratic_eigenvalue_derivative_matches_the_closed_form() {
    // T(λ; b) = λ² I + λ (a I + b L) + L ; mode k: λ² + (a + bμ)λ + μ = 0
    let n = 30;
    let (a, b0, k) = (0.3, 0.02, 4usize);
    let mu = laplacian_1d_eigenvalues(n)[k - 1];
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let build = |b: f64| {
        NonlinearEigenProblem::new(vec![
            NepTerm { matrix: l.clone(), f: Box::new(move |z| C64::ONE + z * c(b, 0.0)) },
            NepTerm { matrix: id.clone(), f: Box::new(move |z| z * z + z * c(a, 0.0)) },
        ])
        .unwrap()
    };
    let p = build(b0);
    let dt = NonlinearEigenProblem::new(vec![NepTerm { matrix: l.clone(), f: Box::new(|z| z) }]).unwrap();
    let d = a + b0 * mu;
    let disc = c(d * d - 4.0 * mu, 0.0).sqrt();
    for sign in [1.0, -1.0] {
        let lam = (c(-d, 0.0) + disc.scale(sign)) * 0.5;
        let x: Vec<C64> =
            (1..=n).map(|j| c((k as f64 * std::f64::consts::PI * j as f64 / (n as f64 + 1.0)).sin(), 0.0)).collect();
        assert!(p.relative_residual(lam, &x) < 1e-12);
        let der = nonlinear_eigenvalue_derivative(&p, &dt, lam, &x, &x);
        // dλ/db = −λ μ / (2λ + a + bμ)
        let exact = -lam * mu / (lam * 2.0 + d);
        assert!((der - exact).norm() < 1e-7 * exact.norm(), "{der} vs {exact}");
    }
    let _ = norm(&[1.0]);
}

#[test]
fn tracking_follows_modes_through_an_exact_crossing() {
    // H(p) = Qᵀ diag(1 + p, 3 − p, 5, 7) Q with Q orthogonal: eigenvalues 1+p and 3−p cross at p = 1
    let q = {
        // a fixed rotation mixing all four coordinates (product of Givens rotations)
        let mut m = Dense::<C64>::identity(4);
        for (i, j, th) in [(0usize, 1usize, 0.7f64), (1, 2, 0.4), (2, 3, 1.1), (0, 3, 0.3)] {
            let mut g = Dense::<C64>::identity(4);
            g[(i, i)] = c(th.cos(), 0.0);
            g[(j, j)] = c(th.cos(), 0.0);
            g[(i, j)] = c(th.sin(), 0.0);
            g[(j, i)] = c(-th.sin(), 0.0);
            m = g.matmul(&m);
        }
        m
    };
    let h_of = |p: f64| {
        let dd = Dense::from_fn(4, 4, |i, j| if i == j { c([1.0 + p, 3.0 - p, 5.0, 7.0][i], 0.0) } else { C64::ZERO });
        let m = q.adjoint().matmul(&dd).matmul(&q);
        let trip: Vec<_> = (0..4).flat_map(|i| (0..4).map(move |j| (i, j))).map(|(i, j)| (i, j, m[(i, j)])).collect();
        Csr::from_triplets(4, 4, &trip).unwrap()
    };
    let steps: Vec<f64> = (0..21).map(|k| 0.05 + 0.09 * k as f64).filter(|p| (p - 1.0).abs() > 1e-9).collect();
    let sorted_modes = |p: f64| {
        let (vals, vecs) = eigs(&h_of(p), None);
        let mut idx: Vec<usize> = (0..4).collect();
        idx.sort_by(|&a, &b| vals[a].re.total_cmp(&vals[b].re));
        (idx.iter().map(|&i| vals[i]).collect::<Vec<_>>(), idx.iter().map(|&i| vecs[i].clone()).collect::<Vec<_>>())
    };
    let (_, v0) = sorted_modes(steps[0]);
    let mut tracker = Tracker::new(v0);
    let (mut max_dev, mut naive_dev) = (0.0f64, 0.0f64);
    for &p in &steps[1..] {
        let (mut vals, mut vecs) = sorted_modes(p);
        naive_dev = naive_dev.max((vals[0].re - (1.0 + p)).abs());
        let mm = tracker.step(&mut vals, &mut vecs, None);
        assert!(mm.overlap.iter().all(|&o| o > 0.5), "{:?} at p = {p}", mm.overlap);
        max_dev = max_dev.max((vals[0].re - (1.0 + p)).abs()).max((vals[1].re - (3.0 - p)).abs());
    }
    assert!(max_dev < 1e-9, "tracked branches deviate from the linear ones by {max_dev}");
    assert!(naive_dev > 0.5, "sorting should lose the branch after the crossing (deviation {naive_dev})");
}

#[test]
fn matching_is_optimal_and_flags_ambiguity() {
    // identical vectors in swapped order: perfect matching, overlap 1
    let e = |k: usize| -> Vec<C64> { (0..3).map(|i| if i == k { C64::ONE } else { C64::ZERO }).collect() };
    let old = vec![e(0), e(1), e(2)];
    let new = vec![e(2), e(0), e(1)];
    let mm = match_modes(&old, &new, None, 0.5, 2.0);
    assert_eq!(mm.assignment, vec![1, 2, 0]);
    assert!(mm.overlap.iter().all(|&o| (o - 1.0).abs() < 1e-14) && mm.uncertain.iter().all(|&u| !u));
    // a 45° mixture is ambiguous
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let mix = vec![vec![c(s, 0.0), c(s, 0.0), C64::ZERO], vec![c(s, 0.0), c(-s, 0.0), C64::ZERO], e(2)];
    let mm = match_modes(&old, &mix, None, 0.6, 1.5);
    assert!(mm.uncertain[0] && mm.uncertain[1] && !mm.uncertain[2], "{:?} {:?}", mm.overlap, mm.uncertain);
}
