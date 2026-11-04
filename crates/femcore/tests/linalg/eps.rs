#![allow(clippy::needless_range_loop)] // index loops mirror the mathematics in these tests
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, laplacian_2d, mass_1d};
use femcore::linalg::eps::{EpsOptions, EpsReason, Which, krylov_schur};
use femcore::linalg::{C64, Csr, LinearOperator, Scalar};

/// Every returned pair satisfies the contract: unit vectors with true residual
/// `‖Ax − λx‖ ≤ slack·tol·|λ|`, and (for Hermitian problems) mutually orthogonal.
fn check_pairs<T: Scalar>(
    a: &Csr<T>,
    b: Option<&Csr<T>>,
    values: &[C64],
    vecs: &femcore::linalg::basis::Basis<T>,
    tol: f64,
    slack: f64,
    orthogonal: bool,
) {
    for i in 0..values.len() {
        let x = vecs.col(i);
        if b.is_none() {
            assert!((norm(x) - 1.0).abs() < 1e-10, "vector {i} not unit: {}", norm(x));
        }
        let r = true_residual(a, b, values[i], x);
        assert!(r <= slack * tol * values[i].norm(), "pair {i}: residual {r} vs tol·|λ| = {}", tol * values[i].norm());
        if orthogonal {
            for j in 0..i {
                let d = x.iter().zip(vecs.col(j)).fold(T::ZERO, |s, (p, q)| s + q.conj() * *p);
                assert!(d.modulus() < 1e-7, "vectors {i},{j} not orthogonal: {}", d.modulus());
            }
        }
    }
}

#[test]
fn laplacian_1d_real_largest_and_smallest() {
    let n = 300;
    let a = laplacian_1d(n);
    let exact = laplacian_1d_eigenvalues(n);
    for (which, nev, targets) in [
        (Which::LargestMagnitude, 5usize, (0..5).map(|i| exact[n - 1 - i]).collect::<Vec<_>>()),
        (Which::SmallestReal, 4, exact[..4].to_vec()),
    ] {
        let opts = EpsOptions { which, tol: 1e-10, max_it: Some(5000), ..EpsOptions::hermitian(nev) };
        let res = krylov_schur(&a, None, None, &opts).unwrap();
        assert_eq!(res.reason, EpsReason::Converged, "{which:?}");
        assert!(res.nconv >= nev);
        for i in 0..nev {
            assert!((res.values[i].re - targets[i]).abs() < 1e-9, "{which:?} {i}: {} vs {}", res.values[i], targets[i]);
            assert_eq!(res.values[i].im, 0.0);
        }
        check_pairs(&a, None, &res.values[..nev], &res.vectors, 1e-10, 5.0, true);
    }
}

/// A complex Hermitian tridiagonal matrix `-e^{iφ}` below / `-e^{-iφ}` above the diagonal is a
/// gauge transform of the 1-D Laplacian: same spectrum, genuinely complex vectors.
fn gauge_laplacian(n: usize, phi: f64) -> Csr<C64> {
    let ph = C64::new(phi.cos(), phi.sin());
    let mut t = Vec::new();
    for i in 0..n {
        if i > 0 {
            t.push((i, i - 1, -ph));
        }
        t.push((i, i, C64::new(2.0, 0.0)));
        if i + 1 < n {
            t.push((i, i + 1, -ph.conj()));
        }
    }
    Csr::from_triplets(n, n, &t).unwrap()
}

#[test]
fn complex_hermitian_gauge_laplacian() {
    let n = 250;
    let a = gauge_laplacian(n, 0.7);
    assert!(a.is_hermitian(1e-15) && !a.is_symmetric(1e-3));
    let exact = laplacian_1d_eigenvalues(n);
    let opts = EpsOptions { tol: 1e-10, max_it: Some(5000), ..EpsOptions::hermitian(4) };
    let res = krylov_schur(&a, None, None, &opts).unwrap();
    assert_eq!(res.reason, EpsReason::Converged);
    for i in 0..4 {
        assert!((res.values[i].re - exact[n - 1 - i]).abs() < 1e-9);
    }
    check_pairs(&a, None, &res.values[..4], &res.vectors, 1e-10, 5.0, true);
}

fn assert_valid_top_spectrum(got: &[f64], exact: &[f64], tol: f64) {
    let mut avail: Vec<f64> = exact.to_vec();
    for &g in got {
        let k = avail
            .iter()
            .position(|e| (e - g).abs() < tol)
            .unwrap_or_else(|| panic!("{g} is not an eigenvalue (or exceeds its multiplicity)"));
        avail.remove(k);
    }
    let lowest = got.iter().cloned().fold(f64::INFINITY, f64::min);
    for &e in exact.iter().filter(|&&e| e > lowest + 100.0 * tol) {
        assert!(got.iter().any(|g| (g - e).abs() < tol), "missed distinct eigenvalue {e} (lowest returned {lowest})");
    }
}

#[test]
fn degenerate_spectrum_single_vector_is_never_wrong() {
    // A square grid has doubly degenerate eigenvalues; a single-vector Krylov method may return
    // one copy only, but must never return a wrong value or skip a distinct one.
    let (nx, ny) = (24, 24);
    let a = laplacian_2d(nx, ny);
    let exact = laplacian_2d_eigenvalues(nx, ny);
    let opts = EpsOptions { tol: 1e-10, ..EpsOptions::hermitian(6) };
    let res = krylov_schur(&a, None, None, &opts).unwrap();
    assert_eq!(res.reason, EpsReason::Converged);
    let got: Vec<f64> = res.values.iter().map(|v| v.re).collect();
    assert_valid_top_spectrum(&got, &exact, 1e-8);
}

#[test]
fn block_variant_finds_every_copy_of_degenerate_eigenvalues() {
    let (nx, ny) = (24, 24);
    let a = laplacian_2d(nx, ny);
    let exact = laplacian_2d_eigenvalues(nx, ny);
    let n = nx * ny;
    // The six largest including multiplicity: pairs (λ_ij, λ_ji) are doubly degenerate.
    let want: Vec<f64> = (0..6).map(|i| exact[n - 1 - i]).collect();
    for b in [2usize, 3] {
        let opts = EpsOptions { tol: 1e-10, block_size: b, ..EpsOptions::hermitian(6) };
        let res = krylov_schur(&a, None, None, &opts).unwrap();
        assert_eq!(res.reason, EpsReason::Converged, "b={b}");
        for i in 0..6 {
            assert!((res.values[i].re - want[i]).abs() < 1e-8, "b={b} i={i}: {} vs {}", res.values[i].re, want[i]);
        }
        check_pairs(&a, None, &res.values[..6], &res.vectors, 1e-10, 5.0, true);
    }
}

#[test]
fn generalised_hermitian_with_mass_weighted_inner_product() {
    // A x = λ M x with the exact FE spectrum λ_k = 6(1−cosθ)/(2+cosθ), θ = kπ/(n+1).
    let n = 80;
    let (a, m) = (laplacian_1d(n), mass_1d(n, 1.0));
    let exact: Vec<f64> = (1..=n)
        .map(|k| {
            let c = (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos();
            6.0 * (1.0 - c) / (2.0 + c)
        })
        .collect();
    let op = DenseOp { n, a: dense_b_inv_a(&a, &m) };
    let w: &dyn LinearOperator<f64> = &m;
    let opts = EpsOptions { which: Which::LargestMagnitude, tol: 1e-11, ..EpsOptions::hermitian(4) };
    let res = krylov_schur(&op, Some(w), None, &opts).unwrap();
    assert_eq!(res.reason, EpsReason::Converged);
    for i in 0..4 {
        assert!((res.values[i].re - exact[n - 1 - i]).abs() < 1e-8, "{i}");
    }
    check_pairs(&a, Some(&m), &res.values[..4], &res.vectors, 1e-11, 20.0, false);
    // Hermitian vectors are returned M-orthonormal (mass-normalised modes), not 2-norm unit.
    let mut mx = vec![0.0; n];
    for i in 0..4 {
        m.spmv(res.vector(i), &mut mx);
        for j in 0..=i {
            let d: f64 = mx.iter().zip(res.vector(j)).map(|(p, q)| p * q).sum();
            assert!((d - if i == j { 1.0 } else { 0.0 }).abs() < 1e-8, "M-orthonormality {i},{j}: {d}");
        }
    }
    assert!(res.b_products > 0);
}

/// `tridiag(-1, 2 + iγ, -1)` is complex symmetric (not Hermitian) with spectrum
/// `λ_k(L) + iγ`, so its eigenvalues are known exactly.
fn complex_symmetric_shifted(n: usize, gamma: f64) -> Csr<C64> {
    Csr::axpby(C64::ONE, &laplacian_1d(n).to_complex(), C64::new(0.0, gamma), &Csr::identity(n))
}

#[test]
fn general_kind_complex_symmetric() {
    let n = 200;
    let gamma = 0.35;
    let a = complex_symmetric_shifted(n, gamma);
    assert!(a.is_symmetric(0.0) && !a.is_hermitian(1e-3));
    let exact = laplacian_1d_eigenvalues(n);
    let opts = EpsOptions { tol: 1e-10, max_it: Some(2000), ..EpsOptions::general(4) };
    let res = krylov_schur(&a, None, None, &opts).unwrap();
    assert_eq!(res.reason, EpsReason::Converged);
    for i in 0..4 {
        let want = C64::new(exact[n - 1 - i], gamma);
        assert!((res.values[i] - want).norm() < 1e-8, "{i}: {} vs {want}", res.values[i]);
    }
    check_pairs(&a, None, &res.values[..4], &res.vectors, 1e-10, 10.0, false);
}

#[test]
fn general_kind_nonnormal_with_known_real_spectrum() {
    // tridiag(-1-c, 2, -1+c) has eigenvalues 2 - 2·sqrt(1-c²)·cos(kπ/(n+1)). Its eigenvalues are
    // exponentially ill-conditioned in n (κ ≈ ((1+c)/(1-c))^{n/2}); n = 80, c = 0.15 keeps
    // κ ≈ 2e5, so the exact spectrum is still a usable reference at 1e-9.
    let (n, c) = (80usize, 0.15f64);
    let mut t = Vec::new();
    for i in 0..n {
        if i > 0 {
            t.push((i, i - 1, C64::new(-1.0 - c, 0.0)));
        }
        t.push((i, i, C64::new(2.0, 0.0)));
        if i + 1 < n {
            t.push((i, i + 1, C64::new(-1.0 + c, 0.0)));
        }
    }
    let a = Csr::from_triplets(n, n, &t).unwrap();
    let s = (1.0 - c * c).sqrt();
    let exact: Vec<f64> =
        (1..=n).map(|k| 2.0 - 2.0 * s * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos()).collect();
    let opts = EpsOptions { which: Which::LargestReal, tol: 1e-9, max_it: Some(3000), ..EpsOptions::general(3) };
    let res = krylov_schur(&a, None, None, &opts).unwrap();
    assert_eq!(res.reason, EpsReason::Converged);
    for i in 0..3 {
        assert!(
            (res.values[i].re - exact[n - 1 - i]).abs() < 1e-6 && res.values[i].im.abs() < 1e-6,
            "{i}: {}",
            res.values[i]
        );
    }
    check_pairs(&a, None, &res.values[..3], &res.vectors, 1e-9, 50.0, false);
}

#[test]
fn general_kind_block_variant_and_history() {
    let n = 150;
    let a = complex_symmetric_shifted(n, 0.2);
    let exact = laplacian_1d_eigenvalues(n);
    let opts =
        EpsOptions { tol: 1e-9, block_size: 2, record_history: true, max_it: Some(3000), ..EpsOptions::general(4) };
    let res = krylov_schur(&a, None, None, &opts).unwrap();
    assert_eq!(res.reason, EpsReason::Converged);
    assert_eq!(res.history.len(), res.iterations);
    assert!(res.history.windows(2).all(|w| w[0].nconv <= w[1].nconv));
    for i in 0..4 {
        assert!((res.values[i] - C64::new(exact[n - 1 - i], 0.2)).norm() < 1e-7);
    }
    assert!(res.timings.total >= res.timings.operator);
}

#[test]
fn real_scalar_rejects_general_kind() {
    let a = laplacian_1d(20);
    assert!(krylov_schur(&a, None, None, &EpsOptions::general(2)).is_err());
    assert!(krylov_schur(&a, None, None, &EpsOptions::hermitian(0)).is_err());
    assert!(krylov_schur(&laplacian_1d(5), None, None, &EpsOptions::hermitian(4)).is_err());
}

#[test]
fn breakdown_in_invariant_subspace_continues() {
    // Diagonal matrix with 3 distinct eigenvalues: the Krylov space of any start vector has
    // dimension 3, so the iteration hits breakdowns and must still return the right values.
    let n = 40;
    let t: Vec<_> = (0..n).map(|i| (i, i, [1.0, 2.0, 3.0][i % 3])).collect();
    let a = Csr::from_triplets(n, n, &t).unwrap();
    let opts = EpsOptions { tol: 1e-10, ..EpsOptions::hermitian(3) };
    let res = krylov_schur(&a, None, None, &opts).unwrap();
    assert_eq!(res.reason, EpsReason::Converged);
    let got: Vec<f64> = res.values.iter().take(3).map(|v| v.re).collect();
    for (g, w) in got.iter().zip([3.0, 3.0, 3.0]) {
        assert!((g - w).abs() < 1e-8, "{got:?}");
    }
    assert!(res.breakdowns > 0);
}
