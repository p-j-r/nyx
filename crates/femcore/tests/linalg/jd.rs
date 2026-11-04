#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, laplacian_2d, mass_1d};
use femcore::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolvePreconditioner, ShiftedSolver};
use femcore::linalg::jd::{JdOptions, jacobi_davidson};
use femcore::linalg::linop::Identity;
use femcore::linalg::{C64, Csr};

fn fe_spectrum(n: usize) -> Vec<f64> {
    (1..=n)
        .map(|k| {
            let c = (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos();
            6.0 * (1.0 - c) / (2.0 + c)
        })
        .collect()
}

fn nearest_set(exact: &[C64], tau: C64, k: usize) -> Vec<C64> {
    let mut v = exact.to_vec();
    v.sort_by(|a, b| (*a - tau).norm().total_cmp(&(*b - tau).norm()));
    v.truncate(k);
    v.sort_by(|a, b| (a.re, a.im).partial_cmp(&(b.re, b.im)).unwrap());
    v
}

/// JD converges to *some* eigenvalues near the target, not necessarily the exact nearest set:
/// each returned value must be a distinct true eigenvalue, among the `slack` nearest.
fn assert_near_true_eigenvalues(got: &[C64], exact: &[C64], tau: C64, slack: usize) {
    let mut by_dist = exact.to_vec();
    by_dist.sort_by(|a, b| (*a - tau).norm().total_cmp(&(*b - tau).norm()));
    let mut used = vec![false; by_dist.len()];
    for g in got {
        let (idx, d) =
            by_dist.iter().enumerate().map(|(i, e)| (i, (*g - *e).norm())).min_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
        assert!(d < 1e-7, "{g} is not an eigenvalue (closest off by {d:e})");
        assert!(idx < slack, "{g} has rank {idx} by distance to the target (allowed {slack})");
        // several copies of a multiple eigenvalue are fine, but each copy only once
        let copy = (0..by_dist.len())
            .find(|&i| !used[i] && (by_dist[i] - *g).norm() < 1e-7)
            .expect("eigenvalue reported more often than its multiplicity");
        used[copy] = true;
    }
}

fn sorted(v: &[C64]) -> Vec<C64> {
    let mut v = v.to_vec();
    v.sort_by(|a, b| (a.re, a.im).partial_cmp(&(b.re, b.im)).unwrap());
    v
}

fn shifted_pc(a: &Csr<C64>, b: Option<&Csr<C64>>, shift: C64) -> ShiftedFactorization<C64> {
    let mut f = ShiftedFactorization::new(a, b, Definiteness::General).unwrap();
    f.set_shift(shift).unwrap();
    f
}

#[test]
fn hermitian_generalised_interior_modes_with_an_exact_shifted_preconditioner() {
    let n = 1500;
    let (a, b) = (laplacian_1d(n).to_complex(), mass_1d(n, 1.0).to_complex());
    let exact: Vec<C64> = fe_spectrum(n).into_iter().map(|x| C64::new(x, 0.0)).collect();
    let tau = C64::new(1.2345, 0.0);
    let f = shifted_pc(&a, Some(&b), tau);
    let pc = ShiftedSolvePreconditioner::new(&f);
    let r = jacobi_davidson(&a, Some(&b), &pc, None, &JdOptions::new(4, tau, true)).unwrap();
    assert!(r.converged, "{} iterations, residuals {:?}", r.iterations, r.history.last());
    assert_near_true_eigenvalues(&r.values, &exact, tau, 12);
    assert!(r.residuals.iter().all(|&x| x < 1e-7), "{:?}", r.residuals);
    // B-orthonormal eigenvectors.
    let mut bx = vec![C64::ZERO; n];
    for i in 0..4 {
        b.spmv(r.vectors.col(i), &mut bx);
        for j in 0..=i {
            let d: C64 = r.vectors.col(j).iter().zip(&bx).map(|(p, q)| p.conj() * *q).sum();
            assert!((d - if i == j { C64::ONE } else { C64::ZERO }).norm() < 1e-5, "B-orthonormality {i},{j}: {d}");
        }
    }
}

#[test]
fn inexact_preconditioners_still_converge() {
    let n = 1500;
    let (a, b) = (laplacian_1d(n).to_complex(), mass_1d(n, 1.0).to_complex());
    let exact: Vec<C64> = fe_spectrum(n).into_iter().map(|x| C64::new(x, 0.0)).collect();
    let tau = C64::new(1.2345, 0.0);
    let mut iters = Vec::new();
    // (A − σB)⁻¹ at a displaced real shift, and at complex-shifted ones: increasingly crude
    // approximations of (A − τB)⁻¹.
    for shift in [tau + C64::new(0.002, 0.0), tau + C64::new(0.0, 0.05), tau + C64::new(0.0, 0.3)] {
        let f = shifted_pc(&a, Some(&b), shift);
        let pc = ShiftedSolvePreconditioner::new(&f);
        let r = jacobi_davidson(&a, Some(&b), &pc, None, &JdOptions::new(4, tau, true)).unwrap();
        assert!(r.converged, "shift {shift}: not converged after {} its", r.iterations);
        assert_near_true_eigenvalues(&r.values, &exact, tau, 12);
        iters.push(r.iterations);
    }
    // The cruder the preconditioner, the more outer iterations.
    assert!(iters[0] <= iters[2] + 2, "{iters:?}");
}

#[test]
fn standard_non_hermitian_complex_symmetric_uses_schur_deflation() {
    // A = L + iγI: eigenvalues λ_k(L) + iγ, eigenvectors real; target between two of them.
    let (n, gamma) = (1200usize, 0.3);
    let a = Csr::axpby(C64::ONE, &laplacian_1d(n).to_complex(), C64::new(0.0, gamma), &Csr::identity(n));
    let ev = laplacian_1d_eigenvalues(n);
    let exact: Vec<C64> = ev.iter().map(|&e| C64::new(e, gamma)).collect();
    let tau = C64::new(0.8, gamma + 0.002);
    let f = shifted_pc(&a, None, tau + C64::new(0.0, 0.01));
    let pc = ShiftedSolvePreconditioner::new(&f);
    let r = jacobi_davidson(&a, None, &pc, None, &JdOptions::new(4, tau, false)).unwrap();
    assert!(r.converged, "{} its", r.iterations);
    for (g, w) in sorted(&r.values).iter().zip(&nearest_set(&exact, tau, 4)) {
        assert!((*g - *w).norm() < 1e-7, "{g} vs {w}");
    }
    assert!(r.residuals.iter().all(|&x| x < 1e-6), "{:?}", r.residuals);
    for i in 0..4 {
        assert!((norm(r.vectors.col(i)) - 1.0).abs() < 1e-10);
    }
}

#[test]
fn degenerate_eigenvalues_are_found_by_deflation() {
    // A square grid: every eigenvalue near 3 is (at least) double. Deflation of the converged
    // vector lets the second copy emerge.
    let (nx, ny) = (24usize, 24usize);
    let a = laplacian_2d(nx, ny).to_complex();
    let exact: Vec<C64> = laplacian_2d_eigenvalues(nx, ny).into_iter().map(|x| C64::new(x, 0.0)).collect();
    let tau = C64::new(3.0013, 0.0);
    let f = shifted_pc(&a, None, tau);
    let pc = ShiftedSolvePreconditioner::new(&f);
    let r = jacobi_davidson(&a, None, &pc, None, &JdOptions::new(4, tau, true)).unwrap();
    assert!(r.converged);
    assert_near_true_eigenvalues(&r.values, &exact, tau, 16);
    // deflation must deliver both copies of the double eigenvalue, not the same one twice
    let at3 = r.values.iter().filter(|v| (v.re - 3.0).abs() < 1e-7).count();
    assert!(at3 == 2 || at3 == 0, "{:?}", r.values);
}

#[test]
fn invalid_and_unsupported_inputs() {
    let n = 100;
    let a = laplacian_1d(n).to_complex();
    let b = mass_1d(n, 1.0).to_complex();
    let tau = C64::new(0.5, 0.0);
    // non-Hermitian pencil
    assert!(jacobi_davidson(&a, Some(&b), &Identity, None, &JdOptions::new(2, tau, false)).is_err());
    // bad dimensions
    assert!(jacobi_davidson(&a, None, &Identity, None, &JdOptions::new(0, tau, true)).is_err());
    assert!(
        jacobi_davidson(&a, None, &Identity, None, &JdOptions { m_max: 4, ..JdOptions::new(2, tau, true) }).is_err()
    );
}
