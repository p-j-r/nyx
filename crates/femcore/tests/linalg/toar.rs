//! TOAR for the quadratic eigenproblem: closed form for modally damped vibration, a general
//! (non-modal) problem against Arnoldi on the explicit linearisation, and the memory saving.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, mass_1d};
use femcore::linalg::dense::Dense;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::faer_glue::generalized_eig_c64;
use femcore::linalg::problem::{EigenProblem, KrylovSchurShiftInvert, ShiftInvertEigenSolver, Symmetry};
use femcore::linalg::toar::{QuadraticProblem, ToarOptions, toar_shift_invert};
use femcore::linalg::{C64, Csr};

fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

fn modal_problem(n: usize, a: f64, b: f64) -> (QuadraticProblem, Vec<C64>) {
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let cm = Csr::axpby(c(a, 0.0), &id, c(b, 0.0), &l);
    let mut all = Vec::new();
    for mu in laplacian_1d_eigenvalues(n) {
        let d = a + b * mu;
        let disc = c(d * d - 4.0 * mu, 0.0).sqrt();
        all.push((c(-d, 0.0) + disc) * 0.5);
        all.push((c(-d, 0.0) - disc) * 0.5);
    }
    (QuadraticProblem { m: id, c: cm, k: l }, all)
}

fn nearest(all: &[C64], sigma: C64, k: usize) -> Vec<C64> {
    let mut v = all.to_vec();
    v.sort_by(|a, b| (*a - sigma).norm().total_cmp(&(*b - sigma).norm()));
    v.truncate(k);
    v
}

fn assert_sets(got: &[C64], want: &[C64], tol: f64) {
    for g in got {
        let d = want.iter().map(|w| (*g - *w).norm()).fold(f64::INFINITY, f64::min);
        assert!(d < tol, "{g} is {d:e} from the closest reference value");
    }
    for w in want {
        let d = got.iter().map(|g| (*g - *w).norm()).fold(f64::INFINITY, f64::min);
        assert!(d < tol, "reference {w} is {d:e} from the closest computed value");
    }
}

#[test]
fn toar_matches_the_closed_form_spectrum_of_modal_damping() {
    let (p, all) = modal_problem(400, 0.3, 0.02);
    let sigma = c(-0.2, 1.0);
    let r = toar_shift_invert(&p, &ToarOptions::new(8, 60, sigma)).unwrap();
    assert!(r.converged, "estimates {:?}", r.estimates);
    assert_sets(&r.values, &nearest(&all, sigma, 8), 1e-8);
    assert!(r.residuals.iter().all(|&x| x < 1e-9), "{:?}", r.residuals);
    for i in 0..8 {
        assert!((norm(r.vectors.col(i)) - 1.0).abs() < 1e-12);
    }
}

#[test]
fn toar_restarts_in_compressed_form_and_halves_the_basis_memory() {
    let (p, all) = modal_problem(2000, 0.3, 0.02);
    let r = toar_shift_invert(&p, &ToarOptions::new(6, 50, c(-0.2, 1.0))).unwrap();
    // this run needs Krylov–Schur restarts (in the compressed form): the result must still be right
    assert!(r.restarts > 0, "expected restarts, got none");
    assert!(r.converged, "estimates {:?}", r.estimates);
    assert_sets(&r.values, &nearest(&all, c(-0.2, 1.0), 6), 1e-8);
    assert!(r.residuals.iter().all(|&x| x < 1e-9), "{:?}", r.residuals);
    let ratio = r.basis_scalars as f64 / r.basis_scalars_linearised as f64;
    eprintln!(
        "TOAR basis: {} scalars vs {} for the linearisation ({ratio:.3}); q = {} for {} steps",
        r.basis_scalars, r.basis_scalars_linearised, r.q_columns, r.steps
    );
    assert!(ratio < 0.55, "ratio {ratio}");
    eprintln!("restarts {}, operator applications {}", r.restarts, r.operator_applications);
}

#[test]
fn toar_agrees_with_arnoldi_on_the_explicit_linearisation_for_a_general_problem() {
    // M = FE mass, K = Laplacian, C = a·(sparse non-modal damping) + gyroscopic-like skew part.
    let n = 300;
    let m = mass_1d(n, 1.0 / n as f64).to_complex();
    let k = laplacian_1d(n).map(|v| v * n as f64).to_complex();
    let mut ct = Vec::new();
    for i in 0..n {
        let d = 0.05 + 0.4 * ((i as f64 / n as f64) - 0.5).powi(2); // position dependent damping
        ct.push((i, i, c(d, 0.0)));
        if i + 1 < n {
            ct.push((i, i + 1, c(0.0, 0.01)));
            ct.push((i + 1, i, c(0.0, -0.01)));
        }
    }
    let cm = Csr::from_triplets(n, n, &ct).unwrap();
    let p = QuadraticProblem { m, c: cm, k };
    let sigma = c(-0.3, 40.0);
    let toar = toar_shift_invert(&p, &ToarOptions { tol: 1e-10, ..ToarOptions::new(6, 80, sigma) }).unwrap();
    assert!(toar.converged, "estimates {:?}", toar.estimates);
    // reference: Krylov–Schur shift-and-invert on the 2n pencil
    let (a, b) = p.linearisation();
    let prob = EigenProblem { a, b: Some(b), symmetry: Symmetry::General };
    let mut ks = KrylovSchurShiftInvert::new(&prob, Definiteness::General).unwrap();
    ks.opts.tol = 1e-12;
    ks.opts.ncv = Some(60);
    let r = ks.solve_near(sigma, 6).unwrap();
    assert!(r.converged);
    // every TOAR pair is one of the converged Krylov–Schur values, and they are the nearest ones
    for g in &toar.values {
        let d = r.values.iter().map(|w| (*g - *w).norm()).fold(f64::INFINITY, f64::min);
        assert!(d < 1e-8, "{g} is {d:e} from the closest Krylov–Schur value");
    }
    assert_sets(&toar.values, &nearest(&r.values, sigma, 6), 1e-8);
    assert!(toar.residuals.iter().all(|&x| x < 1e-8), "{:?}", toar.residuals);
}

#[test]
fn a_small_dense_check_of_the_operator_formula() {
    // n = 6 random-ish problem: all 12 eigenvalues from the companion pencil appear among the Ritz
    // values when ncv covers the whole space.
    let n = 6;
    let t = |f: &dyn Fn(usize, usize) -> C64| {
        let mut tr = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let v = f(i, j);
                if v.norm() > 0.0 {
                    tr.push((i, j, v));
                }
            }
        }
        Csr::from_triplets(n, n, &tr).unwrap()
    };
    let m = t(&|i, j| {
        if i == j {
            c(1.0 + 0.1 * i as f64, 0.0)
        } else if i.abs_diff(j) == 1 {
            c(0.1, 0.0)
        } else {
            C64::new(0.0, 0.0)
        }
    });
    let cc = t(&|i, j| {
        if i == j {
            c(0.2, 0.05 * i as f64)
        } else if i.abs_diff(j) == 1 {
            c(-0.05, 0.02)
        } else {
            C64::new(0.0, 0.0)
        }
    });
    let k = t(&|i, j| {
        if i == j {
            c(2.0 + i as f64, 0.0)
        } else if i.abs_diff(j) == 1 {
            c(-1.0, 0.0)
        } else {
            C64::new(0.0, 0.0)
        }
    });
    let p = QuadraticProblem { m, c: cc, k };
    let (a, b) = p.linearisation();
    let (ad, bd) = (a.to_dense(), b.to_dense());
    let nn = 2 * n;
    let (ev, _) = generalized_eig_c64(
        &Dense::from_fn(nn, nn, |i, j| ad[i * nn + j]),
        &Dense::from_fn(nn, nn, |i, j| bd[i * nn + j]),
    )
    .unwrap();
    // ncv = 12 = 2n: the Krylov space is the whole linearised space, all 12 eigenvalues are exact
    let r = toar_shift_invert(&p, &ToarOptions { tol: 1e-9, ..ToarOptions::new(12, 12, c(0.1, 0.3)) }).unwrap();
    assert!(r.converged, "{:?}", r.estimates);
    assert_sets(&r.values, &ev, 1e-8);
    assert!(r.residuals.iter().all(|&x| x < 1e-8), "{:?}", r.residuals);
}

#[test]
fn invalid_inputs() {
    let (p, _) = modal_problem(20, 0.3, 0.02);
    assert!(toar_shift_invert(&p, &ToarOptions::new(0, 10, c(0.0, 1.0))).is_err());
    assert!(toar_shift_invert(&p, &ToarOptions::new(5, 4, c(0.0, 1.0))).is_err());
    assert!(toar_shift_invert(&p, &ToarOptions::new(5, 41, c(0.0, 1.0))).is_err());
}
