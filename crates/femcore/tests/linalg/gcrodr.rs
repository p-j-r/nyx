// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, laplacian_2d};
use femcore::linalg::ksp::{GcroDr, GcroDrOptions, Gmres, GmresOptions, PcSide, Tolerances};
use femcore::linalg::linop::Identity;
use femcore::linalg::pc::Jacobi;
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Csr, Scalar};

fn rel_residual<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::ZERO; b.len()];
    a.spmv(x, &mut ax);
    ax.iter().zip(b).map(|(p, q)| (*q - *p).modulus_sq()).sum::<f64>().sqrt() / norm(b)
}

fn tol(rtol: f64, max_it: usize) -> Tolerances {
    Tolerances { rtol, max_it, ..Default::default() }
}

/// 1-D convection–diffusion: small eigenvalues of the Laplacian are what stalls restarted GMRES.
fn conv_diff_1d(n: usize, wind: f64) -> Csr<f64> {
    let h = 1.0 / (n as f64 + 1.0);
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 2.0));
        if i > 0 {
            t.push((i, i - 1, -1.0 - wind * h));
        }
        if i + 1 < n {
            t.push((i, i + 1, -1.0 + wind * h));
        }
    }
    Csr::from_triplets(n, n, &t).unwrap()
}

#[test]
fn without_recycling_it_is_right_preconditioned_gmres() {
    let a = conv_diff_1d(200, 20.0);
    let b = random_vector::<f64>(200, 1);
    let jac = Jacobi::new(&a);
    let (mut x1, mut x2) = (vec![0.0; 200], vec![0.0; 200]);
    let t = tol(1e-8, 5000);
    let g = Gmres::new(GmresOptions {
        tol: t,
        restart: 30,
        side: PcSide::Right,
        record_history: true,
        ..Default::default()
    })
    .solve(&a, &jac, &b, &mut x1);
    let r = GcroDr::new(GcroDrOptions { tol: t, m: 30, k: 0, record_history: true, ..Default::default() })
        .solve(&a, &jac, &b, &mut x2);
    assert!(g.reason.is_converged() && r.reason.is_converged());
    assert_eq!(g.iterations, r.iterations);
    for (p, q) in g.history.iter().zip(&r.history) {
        assert!((p - q).abs() <= 1e-4 * p.max(1e-12) + 1e-14, "{p} vs {q}");
    }
}

#[test]
fn recycling_beats_restarted_gmres_on_small_eigenvalues() {
    let n = 400;
    let a = conv_diff_1d(n, 4.0);
    let b = random_vector::<f64>(n, 2);
    let t = tol(1e-8, 20000);
    let mut xg = vec![0.0; n];
    let g = Gmres::new(GmresOptions { tol: t, restart: 30, side: PcSide::Right, ..Default::default() })
        .solve(&a, &Identity, &b, &mut xg);
    let mut x = vec![0.0; n];
    let mut solver = GcroDr::new(GcroDrOptions { tol: t, m: 30, k: 10, ..Default::default() });
    let r = solver.solve(&a, &Identity, &b, &mut x);
    assert!(r.reason.is_converged(), "{:?} after {}", r.reason, r.iterations);
    assert!(rel_residual(&a, &x, &b) < 1e-7, "true residual {}", rel_residual(&a, &x, &b));
    println!(
        "GMRES(30): {} its (conv {}), GCRO-DR(30,10): {} its, {} cycles",
        g.iterations,
        g.reason.is_converged(),
        r.iterations,
        solver.stats().cycles
    );
    assert!(
        !g.reason.is_converged() || r.iterations * 2 < g.iterations,
        "GCRO-DR {} vs GMRES(30) {}",
        r.iterations,
        g.iterations
    );
    let (orth, link) = solver.recycle_defects(&a, &Identity);
    assert!(orth < 1e-9 && link < 1e-8, "orth {orth} link {link}");
    assert!(solver.recycled_dim() > 0);
}

#[test]
fn a_second_right_hand_side_starts_deflated() {
    let n = 400;
    let a = conv_diff_1d(n, 4.0);
    let t = tol(1e-8, 20000);
    let mut solver = GcroDr::new(GcroDrOptions { tol: t, m: 30, k: 10, ..Default::default() });
    let mut its = Vec::new();
    for seed in 0..4 {
        let b = random_vector::<f64>(n, 10 + seed);
        let mut x = vec![0.0; n];
        let info = solver.solve(&a, &Identity, &b, &mut x);
        assert!(info.reason.is_converged());
        assert!(rel_residual(&a, &x, &b) < 1e-7);
        its.push(info.iterations);
    }
    println!("iterations per right-hand side: {its:?}");
    assert!(its[3] * 2 < its[0] || its[1] < its[0], "no gain from recycling: {its:?}");
    assert!(*its.iter().skip(1).min().unwrap() < its[0]);
}

#[test]
fn a_slowly_changing_matrix_reuses_the_space_after_matrix_changed() {
    let n = 300;
    let base = conv_diff_1d(n, 3.0);
    let eye = Csr::<f64>::identity(n);
    let t = tol(1e-8, 20000);
    let b = random_vector::<f64>(n, 3);
    let mut recycler = GcroDr::new(GcroDrOptions { tol: t, m: 30, k: 10, ..Default::default() });
    let (mut with, mut without) = (0usize, 0usize);
    for step in 0..5 {
        let a = Csr::axpby(1.0, &base, 1e-4 * step as f64, &eye);
        let mut x = vec![0.0; n];
        if step > 0 {
            recycler.matrix_changed();
        }
        let info = recycler.solve(&a, &Identity, &b, &mut x);
        assert!(info.reason.is_converged(), "step {step}: {:?}", info.reason);
        assert!(rel_residual(&a, &x, &b) < 1e-7);
        let (orth, link) = recycler.recycle_defects(&a, &Identity);
        assert!(orth < 1e-9 && link < 1e-7, "step {step}: orth {orth} link {link}");
        with += info.iterations;
        let mut fresh = GcroDr::new(GcroDrOptions { tol: t, m: 30, k: 10, ..Default::default() });
        let mut y = vec![0.0; n];
        without += fresh.solve(&a, &Identity, &b, &mut y).iterations;
    }
    println!("sequence of 5 shifted systems: {with} iterations with recycling, {without} without");
    assert!(with < without, "{with} vs {without}");
}

#[test]
fn complex_symmetric_system_with_jacobi() {
    let nx = 30;
    let l = laplacian_2d(nx, nx).to_complex();
    let t: Vec<_> = (0..nx * nx).map(|p| (p, p, C64::new(-1.8, 0.2 + 0.1 * (p % nx) as f64 / nx as f64))).collect();
    let a = Csr::axpby(C64::ONE, &l, C64::ONE, &Csr::from_triplets(nx * nx, nx * nx, &t).unwrap());
    let b = random_vector::<C64>(nx * nx, 5);
    let jac = Jacobi::new(&a);
    let mut x = vec![C64::ZERO; b.len()];
    let mut solver = GcroDr::new(GcroDrOptions { tol: tol(1e-9, 5000), m: 30, k: 8, ..Default::default() });
    let info = solver.solve(&a, &jac, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    assert!(rel_residual(&a, &x, &b) < 1e-7);
    let (orth, link) = solver.recycle_defects(&a, &jac);
    assert!(orth < 1e-9 && link < 1e-8, "orth {orth} link {link}");
}

#[test]
fn spd_laplacian_zero_rhs_and_nonzero_initial_guess() {
    let a = laplacian_1d(100);
    let mut solver = GcroDr::new(GcroDrOptions { tol: tol(1e-10, 5000), m: 20, k: 5, ..Default::default() });
    let mut x = vec![1.0; 100];
    let z = solver.solve(&a, &Identity, &vec![0.0; 100], &mut x);
    assert!(z.reason.is_converged() && z.iterations == 0 && x.iter().all(|v| *v == 0.0));
    let b = random_vector::<f64>(100, 7);
    let exact = {
        let mut y = vec![0.0; 100];
        let mut s = GcroDr::new(GcroDrOptions { tol: tol(1e-13, 5000), m: 40, k: 0, ..Default::default() });
        s.solve(&a, &Identity, &b, &mut y);
        y
    };
    let mut x: Vec<f64> = exact.iter().map(|v| v + 0.1).collect();
    let mut opts = solver.opts.clone();
    opts.initial_guess_nonzero = true;
    solver.opts = opts;
    let info = solver.solve(&a, &Identity, &b, &mut x);
    assert!(info.reason.is_converged());
    let err = x.iter().zip(&exact).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
    assert!(err < 1e-6, "{err}");
}
