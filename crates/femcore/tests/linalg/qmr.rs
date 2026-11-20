// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::laplacian_2d;
use femcore::linalg::ilu::Ilu;
use femcore::linalg::ksp::{Gmres, GmresOptions, Qmr, QmrOptions, Tolerances};
use femcore::linalg::linop::Identity;
use femcore::linalg::pc::Jacobi;
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Csr, Scalar};

fn conv_diff(nx: usize, wind: f64) -> Csr<f64> {
    let h = 1.0 / (nx as f64 + 1.0);
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let p = j * nx + i;
            t.push((p, p, 4.0));
            if i > 0 {
                t.push((p, p - 1, -1.0 - wind * h));
            }
            if i + 1 < nx {
                t.push((p, p + 1, -1.0 + wind * h));
            }
            if j > 0 {
                t.push((p, p - nx, -1.0));
            }
            if j + 1 < nx {
                t.push((p, p + nx, -1.0));
            }
        }
    }
    Csr::from_triplets(nx * nx, nx * nx, &t).unwrap()
}

fn helmholtz(nx: usize, k2: f64, sigma: f64) -> Csr<C64> {
    let l = laplacian_2d(nx, nx).to_complex();
    let t: Vec<_> = (0..nx * nx)
        .map(|p| {
            let (i, j) = (p % nx, p / nx);
            (p, p, C64::new(-k2, sigma * (1.0 + 0.5 * i as f64 / nx as f64 + 0.25 * j as f64 / nx as f64)))
        })
        .collect();
    Csr::axpby(C64::ONE, &l, C64::ONE, &Csr::from_triplets(nx * nx, nx * nx, &t).unwrap())
}

fn rel_residual<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::ZERO; b.len()];
    a.spmv(x, &mut ax);
    ax.iter().zip(b).map(|(p, q)| (*q - *p).modulus_sq()).sum::<f64>().sqrt() / norm(b)
}

fn tol() -> Tolerances {
    Tolerances { rtol: 1e-9, max_it: 3000, ..Default::default() }
}

#[test]
fn qmr_solves_real_nonsymmetric_systems() {
    let a = conv_diff(40, 15.0);
    let at = a.transpose();
    let b = random_vector::<f64>(a.nrows(), 3);
    let jac = Jacobi::new(&a);
    let mut x = vec![0.0; b.len()];
    let info = Qmr::new(QmrOptions { tol: tol(), record_history: true }).solve(&a, &at, &jac, &jac, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    assert!(rel_residual(&a, &x, &b) < 1e-7, "true residual {}", rel_residual(&a, &x, &b));
    assert_eq!(info.matvecs, 2 * info.iterations, "A and Aᵀ once per iteration");
    // comparable work to full GMRES (which needs the fewest iterations)
    let mut xg = vec![0.0; b.len()];
    let g = Gmres::new(GmresOptions { tol: tol(), restart: 400, ..Default::default() }).solve(&a, &jac, &b, &mut xg);
    assert!(g.reason.is_converged());
    assert!(info.iterations < 3 * g.iterations, "QMR {} vs GMRES {}", info.iterations, g.iterations);
}

#[test]
fn qmr_solves_complex_nonsymmetric_systems() {
    // complex Helmholtz plus a convective term: neither symmetric nor Hermitian
    let nx = 24;
    let base = helmholtz(nx, 1.5, 0.4);
    let mut t = Vec::new();
    for p in 0..nx * nx {
        if p % nx + 1 < nx {
            t.push((p, p + 1, C64::new(0.3, 0.1)));
        }
    }
    let a = Csr::axpby(C64::ONE, &base, C64::ONE, &Csr::from_triplets(nx * nx, nx * nx, &t).unwrap());
    assert!(!a.is_symmetric(1e-12));
    let at = a.transpose();
    let b = random_vector::<C64>(nx * nx, 4);
    let jac = Jacobi::new(&a);
    let mut x = vec![C64::ZERO; b.len()];
    let info = Qmr::new(QmrOptions { tol: tol(), record_history: false }).solve(&a, &at, &jac, &jac, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    assert!(rel_residual(&a, &x, &b) < 1e-7);
}

#[test]
fn symmetric_qmr_uses_one_product_and_reproduces_the_general_iterates() {
    let a = helmholtz(30, 2.0, 0.5);
    assert!(a.is_symmetric(0.0));
    let b = random_vector::<C64>(a.nrows(), 8);
    let jac = Jacobi::new(&a);
    let opts = QmrOptions { tol: tol(), record_history: true };
    let (mut xs, mut xg) = (vec![C64::ZERO; b.len()], vec![C64::ZERO; b.len()]);
    let sym = Qmr::new(opts.clone()).solve_symmetric(&a, &jac, &b, &mut xs);
    let gen_ = Qmr::new(opts).solve(&a, &a, &jac, &jac, &b, &mut xg);
    assert!(sym.reason.is_converged() && gen_.reason.is_converged());
    assert_eq!(sym.matvecs, sym.iterations, "one product with A per iteration");
    assert_eq!(gen_.matvecs, 2 * gen_.iterations);
    // same Krylov space and same normalisation of the right Lanczos vectors: same iterates
    assert!(sym.iterations.abs_diff(gen_.iterations) <= 1, "{} vs {}", sym.iterations, gen_.iterations);
    let k = sym.history.len().min(gen_.history.len()).min(40);
    for i in 0..k {
        let (p, q) = (sym.history[i], gen_.history[i]);
        assert!((p - q).abs() <= 1e-8 * q.max(1e-30), "history differs at {i}: {p} vs {q}");
    }
    assert!(rel_residual(&a, &xs, &b) < 1e-7);
    let diff = xs.iter().zip(&xg).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm(&xg);
    assert!(diff < 1e-6, "{diff}");
}

#[test]
fn symmetric_qmr_with_ilu_on_real_spd_matches_cg_scale() {
    // On SPD systems QMR is mathematically close to CG/MINRES; check only that it converges quickly.
    let a = laplacian_2d(40, 40);
    let b = random_vector::<f64>(1600, 2);
    let ilu = Ilu::ilu_k(&a, 0);
    let mut x = vec![0.0; 1600];
    let info = Qmr::new(QmrOptions { tol: tol(), record_history: false }).solve_symmetric(&a, &ilu, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?}", info.reason);
    assert!(rel_residual(&a, &x, &b) < 1e-7);
    assert!(info.iterations < 120, "{}", info.iterations);
}

#[test]
fn zero_rhs_and_identity_preconditioner() {
    let a = laplacian_2d(10, 10);
    let mut x = vec![1.0; 100];
    let info = Qmr::new(QmrOptions::default()).solve_symmetric(&a, &Identity, &vec![0.0; 100], &mut x);
    assert!(info.reason.is_converged() && info.iterations == 0);
    assert!(x.iter().all(|v| *v == 0.0));
}

#[test]
fn history_is_smooth_where_bicgstab_style_peaks_would_appear() {
    // The quasi-residual of QMR is monotone; the true (recursively updated) residual on this problem
    // never grows by more than a modest factor from one step to the next.
    let a = conv_diff(40, 40.0);
    let at = a.transpose();
    let b = random_vector::<f64>(a.nrows(), 9);
    let jac = Jacobi::new(&a);
    let mut x = vec![0.0; b.len()];
    let info = Qmr::new(QmrOptions { tol: tol(), record_history: true }).solve(&a, &at, &jac, &jac, &b, &mut x);
    assert!(info.reason.is_converged());
    let worst = info.history.windows(2).map(|w| w[1] / w[0]).fold(0.0f64, f64::max);
    assert!(worst < 3.0, "largest one-step residual growth {worst}");
}
