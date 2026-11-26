// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::amg::{Amg, AmgOptions, shifted_laplacian_operator};
use femcore::linalg::csr::{laplacian_2d, laplacian_3d};
use femcore::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::linalg::ksp::{Cg, CgKind, CgOptions, Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances};
use femcore::linalg::linop::{Identity, Preconditioner};
use femcore::linalg::rng::random_vector;
use femcore::linalg::vecops;
use femcore::linalg::{C64, Csr, Scalar};

/// `L − k² + iσ·ramp`: complex symmetric Helmholtz-like matrix on an nx × ny grid.
fn helmholtz(nx: usize, ny: usize, k2: f64, sigma: f64) -> Csr<C64> {
    let t: Vec<_> = (0..nx * ny)
        .map(|p| {
            let (i, j) = (p % nx, p / nx);
            let ramp = 1.0 + 0.5 * (i as f64 / nx as f64) + 0.25 * (j as f64 / ny as f64);
            (p, p, C64::new(-k2, sigma * ramp))
        })
        .collect();
    Csr::axpby(
        C64::ONE,
        &laplacian_2d(nx, ny).to_complex(),
        C64::ONE,
        &Csr::from_triplets(nx * ny, nx * ny, &t).unwrap(),
    )
}

fn residual_norm<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::ZERO; b.len()];
    a.spmv(x, &mut ax);
    ax.iter().zip(b).map(|(p, q)| (*q - *p).modulus_sq()).sum::<f64>().sqrt()
}

#[test]
fn hierarchy_coarsens_with_bounded_complexity() {
    for a in [laplacian_2d(64, 64), laplacian_3d(16)] {
        let amg = Amg::new(&a, AmgOptions::default()).unwrap();
        let st = amg.stats();
        assert!(st.level_sizes.windows(2).all(|w| w[1] < w[0]), "{:?}", st.level_sizes);
        assert!(*st.level_sizes.last().unwrap() <= 150);
        assert!(st.operator_complexity < 2.0, "operator complexity {}", st.operator_complexity);
    }
}

#[test]
fn v_cycle_is_a_convergent_stationary_iteration_independent_of_the_grid() {
    let mut rates = Vec::new();
    for n in [32usize, 64, 128] {
        let a = laplacian_2d(n, n);
        let amg = Amg::new(&a, AmgOptions { sweeps: 2, ..Default::default() }).unwrap();
        let b = random_vector::<f64>(n * n, 5);
        let mut x = vec![0.0; n * n];
        let mut r = b.clone();
        let (mut z, r0) = (vec![0.0; n * n], norm(&b));
        for _ in 0..10 {
            amg.apply(&r, &mut z);
            vecops::axpy(1.0, &z, &mut x);
            r = {
                let mut ax = vec![0.0; n * n];
                a.spmv(&x, &mut ax);
                b.iter().zip(&ax).map(|(p, q)| p - q).collect()
            };
        }
        rates.push((norm(&r) / r0).powf(0.1));
    }
    assert!(rates.iter().all(|&r| r < 0.4), "convergence factors {rates:?}");
    assert!(rates[2] < 1.5 * rates[0] + 0.05, "factors should not degrade with the grid: {rates:?}");
}

#[test]
fn cycle_is_symmetric_for_complex_symmetric_matrices() {
    let a = helmholtz(40, 40, -0.5, 0.7); // positive real shift, absorbing
    let amg = Amg::new(&a, AmgOptions::default()).unwrap();
    let n = a.nrows();
    let (x, y) = (random_vector::<C64>(n, 1), random_vector::<C64>(n, 2));
    let (mut mx, mut my) = (vec![C64::ZERO; n], vec![C64::ZERO; n]);
    amg.apply(&x, &mut mx);
    amg.apply(&y, &mut my);
    // Symmetric in the bilinear form: yᵀ M x = xᵀ M y.
    let (l, r) = (vecops::dotu(&y, &mx), vecops::dotu(&x, &my));
    assert!((l - r).norm() < 1e-11 * l.norm(), "{l} vs {r}");
}

#[test]
fn pcg_with_amg_has_grid_independent_iteration_counts() {
    let mut its = Vec::new();
    for n in [32usize, 64, 128] {
        let a = laplacian_2d(n, n);
        let amg = Amg::new(&a, AmgOptions::default()).unwrap();
        let b = random_vector::<f64>(n * n, 7);
        let mut x = vec![0.0; n * n];
        let info = Cg::new(CgOptions { tol: Tolerances { rtol: 1e-8, ..Default::default() }, ..Default::default() })
            .solve(&a, &amg, &b, &mut x);
        assert!(info.reason.is_converged());
        its.push(info.iterations);
    }
    assert!(its.iter().all(|&i| i <= 25), "{its:?}");
    assert!(its[2] <= its[0] + 6, "{its:?}");
}

#[test]
fn cocg_with_amg_solves_a_complex_symmetric_problem_in_few_iterations() {
    let mut its = Vec::new();
    for n in [40usize, 80, 160] {
        let a = helmholtz(n, n, -0.05, 0.3);
        let amg = Amg::new(&a, AmgOptions::default()).unwrap();
        let b = random_vector::<C64>(n * n, 3);
        let mut x = vec![C64::ZERO; n * n];
        let opts = CgOptions {
            tol: Tolerances { rtol: 1e-8, ..Default::default() },
            kind: CgKind::Symmetric,
            ..Default::default()
        };
        let info = Cg::new(opts).solve(&a, &amg, &b, &mut x);
        assert!(info.reason.is_converged(), "n={n}: {:?}", info.reason);
        assert!(residual_norm(&a, &x, &b) < 1e-7 * norm(&b));
        its.push(info.iterations);
    }
    assert!(its.iter().all(|&i| i <= 30), "{its:?}");
    assert!(its[2] <= its[0] + 8, "{its:?}");
}

#[test]
fn hermitian_option_gives_a_working_preconditioner_for_cg() {
    let n = 80;
    let ph = C64::new(0.6f64.cos(), 0.6f64.sin());
    // 2-D complex Hermitian: Laplacian with a phase on the x-links, plus a positive shift.
    let mut t = Vec::new();
    for j in 0..n {
        for i in 0..n {
            let p = i + n * j;
            t.push((p, p, C64::new(4.1, 0.0)));
            if i > 0 {
                t.push((p, p - 1, -ph));
            }
            if i + 1 < n {
                t.push((p, p + 1, -ph.conj()));
            }
            if j > 0 {
                t.push((p, p - n, C64::new(-1.0, 0.0)));
            }
            if j + 1 < n {
                t.push((p, p + n, C64::new(-1.0, 0.0)));
            }
        }
    }
    let h = Csr::from_triplets(n * n, n * n, &t).unwrap();
    assert!(h.is_hermitian(1e-14));
    let amg = Amg::new(&h, AmgOptions { hermitian: true, ..Default::default() }).unwrap();
    let b = random_vector::<C64>(n * n, 4);
    let mut x = vec![C64::ZERO; n * n];
    let info = Cg::new(CgOptions { tol: Tolerances { rtol: 1e-8, ..Default::default() }, ..Default::default() })
        .solve(&h, &amg, &b, &mut x);
    assert!(info.reason.is_converged() && info.iterations <= 30, "{info:?}");
}

/// Exact solve of the shifted operator, the best case any approximation of it can reach.
struct ExactShifted(ShiftedFactorization<C64>);

impl Preconditioner<C64> for ExactShifted {
    fn apply(&self, r: &[C64], z: &mut [C64]) {
        z.copy_from_slice(r);
        self.0.solve_block_in_place(z, 1).unwrap();
    }
}

#[test]
fn shifted_laplacian_amg_is_as_good_as_an_exact_shifted_solve_and_far_better_than_nothing() {
    // kh = 0.4 → ~16 points per wavelength, light damping: the real part of A has many negative
    // eigenvalues, so COCG fails and AMG on A itself is useless. One AMG V-cycle on the shifted
    // Laplacian K − (1+0.5i)k²I works about as well as inverting that operator exactly.
    let k2 = 0.16;
    let beta = C64::new(1.0, 0.5);
    for n in [40usize, 60] {
        let a = helmholtz(n, n, k2, 0.01);
        let shifted = shifted_laplacian_operator(&laplacian_2d(n, n), None, k2, beta);
        let amg = Amg::new(&shifted, AmgOptions::default()).unwrap();
        let mut f = ShiftedFactorization::new(&shifted, None, Definiteness::General).unwrap();
        f.set_shift(C64::ZERO).unwrap();
        let exact = ExactShifted(f);
        let b = random_vector::<C64>(n * n, 11);
        let tol = Tolerances { rtol: 1e-8, max_it: 3000, ..Default::default() };
        let opts = GmresOptions {
            tol,
            restart: 100,
            side: PcSide::Right,
            flexible: true,
            orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
            ..Default::default()
        };
        let run = |pc: &dyn Preconditioner<C64>| {
            let mut x = vec![C64::ZERO; n * n];
            let info = Gmres::new(opts.clone()).solve(&a, pc, &b, &mut x);
            (info, x)
        };
        let (with_amg, x) = run(&amg);
        let (with_exact, _) = run(&exact);
        assert!(with_amg.reason.is_converged() && with_exact.reason.is_converged(), "n={n}");
        assert!(residual_norm(&a, &x, &b) < 1e-7 * norm(&b));
        assert!(
            with_amg.iterations as f64 <= 1.3 * with_exact.iterations as f64 + 5.0,
            "n={n}: AMG {} vs exact shifted solve {}",
            with_amg.iterations,
            with_exact.iterations
        );
        // No preconditioner: GMRES(100) does not get there within five times that many iterations.
        let budget = 5 * with_exact.iterations;
        let mut x0 = vec![C64::ZERO; n * n];
        let plain = Gmres::new(GmresOptions { tol: Tolerances { max_it: budget, ..tol }, ..opts.clone() })
            .solve(&a, &Identity, &b, &mut x0);
        assert!(!plain.reason.is_converged(), "n={n}: plain GMRES converged within {budget} iterations");
    }
}
