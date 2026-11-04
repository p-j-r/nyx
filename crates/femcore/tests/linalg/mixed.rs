//! Mixed-precision LU with refinement: accuracy equal to a double-precision solve, detection of matrices that
//! single precision cannot handle, and use inside a shift-and-invert eigensolve.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, laplacian_2d, laplacian_3d, mass_1d};
use femcore::linalg::direct::{Definiteness, ShiftInvertOperator, ShiftedFactorization, ShiftedSolver};
use femcore::linalg::eps::{EpsOptions, Which, krylov_schur};
use femcore::linalg::mixed::MixedPrecisionLu;
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Csr};
use std::time::Instant;

fn convdiff(nx: usize) -> Csr<f64> {
    let h = 1.0 / (nx as f64 + 1.0);
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let p = j * nx + i;
            t.push((p, p, 4.0));
            if i > 0 {
                t.push((p, p - 1, -1.0 - 20.0 * h));
            }
            if i + 1 < nx {
                t.push((p, p + 1, -1.0 + 20.0 * h));
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

#[test]
fn real_system_is_solved_to_double_precision() {
    let a = convdiff(60);
    let n = a.nrows();
    let mut mixed = MixedPrecisionLu::new(&a, None).unwrap();
    mixed.set_shift(0.0).unwrap();
    let mut exact = ShiftedFactorization::new(&a, None, Definiteness::General).unwrap();
    exact.set_shift(0.0).unwrap();
    let b = random_vector::<f64>(n * 3, 1);
    let (mut xm, mut xe) = (b.clone(), b.clone());
    mixed.solve_block_in_place(&mut xm, 3);
    exact.solve_block_in_place(&mut xe, 3);
    let err = xm.iter().zip(&xe).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&xe);
    let st = mixed.mixed_stats();
    eprintln!(
        "real convection-diffusion: error vs double {err:.2e}, {} refinement steps for 3 solves",
        st.refinement_steps
    );
    assert!(err < 1e-11, "{err}");
    assert_eq!(st.fallback_solves, 0);
    assert!(st.refinement_steps <= 3 * 5, "{}", st.refinement_steps);
}

#[test]
fn complex_shifted_pencil_with_a_shift_and_a_mass_matrix() {
    let n = 40 * 40;
    let a = laplacian_2d(40, 40).to_complex();
    let b: Csr<C64> = Csr::axpby(
        C64::new(1.0, 0.0),
        &Csr::identity(n),
        C64::new(0.0, 0.2),
        &laplacian_2d(40, 40).to_complex().map(|v| v * C64::new(0.01, 0.0)),
    );
    let sigma = C64::new(1.3, 0.4);
    let mut mixed = MixedPrecisionLu::new(&a, Some(&b)).unwrap();
    mixed.set_shift(sigma).unwrap();
    let mut exact = ShiftedFactorization::new(&a, Some(&b), Definiteness::General).unwrap();
    exact.set_shift(sigma).unwrap();
    let rhs = random_vector::<C64>(n, 2);
    let (mut xm, mut xe) = (rhs.clone(), rhs.clone());
    mixed.solve_block_in_place(&mut xm, 1);
    exact.solve_block_in_place(&mut xe, 1);
    let err = xm.iter().zip(&xe).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm(&xe);
    assert!(err < 1e-11, "{err}");
    assert_eq!(mixed.mixed_stats().fallback_solves, 0);
}

#[test]
fn nearly_singular_shift_is_detected_and_handled_by_the_fallback() {
    // shift 1e-9 relative away from an eigenvalue: κ ≈ 1e9 ≫ 1/ε₃₂, single precision is useless
    let n = 200;
    let a = laplacian_1d(n);
    let lam = laplacian_1d_eigenvalues(n)[9];
    let mut mixed = MixedPrecisionLu::new(&a, None).unwrap();
    mixed.set_shift(lam * (1.0 + 1e-9)).unwrap();
    let b = random_vector::<f64>(n, 5);
    let mut x = b.clone();
    mixed.solve_block_in_place(&mut x, 1);
    let st = mixed.mixed_stats();
    assert!(st.fallback_solves > 0 || st.low_factor_failures > 0, "{st:?}");
    // the answer solves the system to the accuracy of the double-precision factorisation
    let mut ax = vec![0.0; n];
    a.spmv(&x, &mut ax);
    for i in 0..n {
        ax[i] -= lam * (1.0 + 1e-9) * x[i];
    }
    let res = ax.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&b);
    // reference: what a double-precision LU achieves on the same matrix (backward stable: residual ~ ε‖A‖‖x‖/‖b‖ ~ 1e-7 at κ = 1e9)
    let sh = Csr::axpby(1.0, &a, -lam * (1.0 + 1e-9), &Csr::identity(n));
    let mut exact = ShiftedFactorization::new(&sh, None, Definiteness::General).unwrap();
    exact.set_shift(0.0).unwrap();
    let mut xe = b.clone();
    exact.solve_block_in_place(&mut xe, 1);
    let mut r2 = vec![0.0; n];
    sh.spmv(&xe, &mut r2);
    let res_exact = r2.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&b);
    assert!(res <= 10.0 * res_exact + 1e-12, "mixed residual {res:e} vs double {res_exact:e}");
}

#[test]
fn shift_and_invert_eigensolve_with_the_mixed_solver_matches_the_exact_one() {
    let n = 600;
    let a = laplacian_1d(n);
    let m = mass_1d(n, 1.0);
    let sigma = 3.3;
    let run = |solver: &dyn femcore::linalg::linop::LinearOperator<f64>| {
        let mut o = EpsOptions::hermitian(4);
        o.which = Which::LargestMagnitude;
        o.tol = 1e-10;
        o.ncv = Some(30);
        let r = krylov_schur(solver, None, None, &o).unwrap();
        let mut v: Vec<f64> = r.values.iter().take(4).map(|t| sigma + 1.0 / t.re).collect();
        v.sort_by(f64::total_cmp);
        v
    };
    let mut mixed = MixedPrecisionLu::new(&a, Some(&m)).unwrap();
    mixed.set_shift(sigma).unwrap();
    let mut exact = ShiftedFactorization::new(&a, Some(&m), Definiteness::General).unwrap();
    exact.set_shift(sigma).unwrap();
    let vm = run(&ShiftInvertOperator::new(&mixed, Some(&m)));
    let ve = run(&ShiftInvertOperator::new(&exact, Some(&m)));
    for (p, q) in vm.iter().zip(&ve) {
        assert!((p - q).abs() < 1e-8 * q.abs().max(1.0), "{p} vs {q}");
    }
    let st = mixed.mixed_stats();
    eprintln!(
        "eigensolve: {} refined solves, {} refinement steps, {} fallbacks",
        st.refined_solves, st.refinement_steps, st.fallback_solves
    );
    // an occasional fallback is legitimate: Krylov vectors that approach an eigenvector make the shifted system more ill-conditioned
    assert!(st.fallback_solves <= 3, "{st:?}");
}

#[test]
#[ignore = "timing study; run with --ignored --nocapture"]
fn timing_3d_real() {
    let a = laplacian_3d(40);
    let shifted = Csr::axpby(1.0, &a, -0.001, &Csr::identity(a.nrows()));
    let t = Instant::now();
    let mut exact = ShiftedFactorization::new(&shifted, None, Definiteness::General).unwrap();
    exact.set_shift(0.0).unwrap();
    let te = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let mut mixed = MixedPrecisionLu::new(&shifted, None).unwrap();
    mixed.set_shift(0.0).unwrap();
    let tm = t.elapsed().as_secs_f64();
    let b = random_vector::<f64>(shifted.nrows(), 3);
    let mut x1 = b.clone();
    let t = Instant::now();
    exact.solve_block_in_place(&mut x1, 1);
    let se = t.elapsed().as_secs_f64();
    let mut x2 = b.clone();
    let t = Instant::now();
    mixed.solve_block_in_place(&mut x2, 1);
    let sm = t.elapsed().as_secs_f64();
    eprintln!(
        "3-D LU n = 64000: factor {te:.2} s (f64) vs {tm:.2} s (f32); solve {se:.3} s vs {sm:.3} s (with refinement, {:?})",
        mixed.mixed_stats()
    );
}
