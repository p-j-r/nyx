//! Nonlinear eigenproblems with closed-form spectra: each problem is a function of the 1-D Laplacian
//! `L` (eigenvalues μ_k, sine eigenvectors), so every mode reduces to a scalar equation whose roots
//! are known: damped vibration (quadratic) and a Drude-type dispersive medium (rational, cubic).
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::contour::Ellipse;
use femcore::linalg::csr::laplacian_1d;
use femcore::linalg::dense::Dense;
use femcore::linalg::faer_glue::generalized_eig_c64;
use femcore::linalg::nep::{
    BeynOptions, NepTerm, NonlinearEigenProblem, RiiOptions, beyn_contour, residual_inverse_iteration,
};
use femcore::linalg::{C64, Csr};

fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

fn roots(poly_low_to_high: &[C64]) -> Vec<C64> {
    // companion matrix eigenvalues of a monic-normalised polynomial
    let d = poly_low_to_high.len() - 1;
    let lead = poly_low_to_high[d];
    let a = Dense::from_fn(d, d, |i, j| {
        if i == 0 {
            -poly_low_to_high[d - 1 - j] / lead
        } else if i == j + 1 {
            C64::ONE
        } else {
            C64::ZERO
        }
    });
    generalized_eig_c64(&a, &Dense::identity(d)).unwrap().0
}

fn sorted(v: &[C64]) -> Vec<C64> {
    let mut v = v.to_vec();
    v.sort_by(|a, b| (a.re, a.im).partial_cmp(&(b.re, b.im)).unwrap());
    v
}

fn assert_same_sets(got: &[C64], want: &[C64], tol: f64) {
    assert_eq!(got.len(), want.len(), "got {got:?}\nwant {want:?}");
    for (g, w) in sorted(got).iter().zip(&sorted(want)) {
        assert!((*g - *w).norm() < tol, "{g} vs {w}");
    }
}

/// Damped vibration `λ² I + λ (a I + b L) + L`.
fn damped(n: usize, a: f64, b: f64) -> (NonlinearEigenProblem, Vec<C64>) {
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let cm = Csr::axpby(c(a, 0.0), &id, c(b, 0.0), &l);
    let p = NonlinearEigenProblem::polynomial(vec![l, cm, id]).unwrap();
    let mut all = Vec::new();
    for mu in laplacian_1d_eigenvalues(n) {
        all.extend(roots(&[c(mu, 0.0), c(a + b * mu, 0.0), C64::ONE]));
    }
    (p, all)
}

/// `L − λ² (1 − ω_p²/(λ² + iγλ)) I = L − λ² I + ω_p² λ/(λ + iγ) I`.
fn drude(n: usize, wp: f64, gamma: f64) -> (NonlinearEigenProblem, Vec<C64>) {
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let p = NonlinearEigenProblem::new(vec![
        NepTerm { matrix: l, f: Box::new(|_| C64::ONE) },
        NepTerm { matrix: id, f: Box::new(move |z| -z * z + c(wp * wp, 0.0) * z / (z + c(0.0, gamma))) },
    ])
    .unwrap();
    let mut all = Vec::new();
    for mu in laplacian_1d_eigenvalues(n) {
        // λ³ + iγλ² − (μ + ω_p²)λ − iγμ = 0
        all.extend(roots(&[c(0.0, -gamma * mu), c(-(mu + wp * wp), 0.0), c(0.0, gamma), C64::ONE]));
    }
    (p, all)
}

fn inside(all: &[C64], region: &Ellipse) -> Vec<C64> {
    all.iter().copied().filter(|z| region.contains(*z)).collect()
}

#[test]
fn beyn_finds_all_quadratic_eigenvalues_in_a_region() {
    let (p, all) = damped(60, 0.3, 0.02);
    let region = Ellipse { center: c(-0.15, 0.4), a: 0.1, b: 0.165 };
    let want = inside(&all, &region);
    assert!(want.len() >= 6, "{}", want.len());
    let r = beyn_contour(&p, &BeynOptions { nodes: 48, ..BeynOptions::new(region) }).unwrap();
    assert_same_sets(&r.values, &want, 1e-7);
    assert!(r.residuals.iter().all(|&x| x < 1e-8), "{:?}", r.residuals);
    assert!(r.rank >= want.len(), "rank {} < {}", r.rank, want.len());
    assert!(!r.saturated);
}

#[test]
fn beyn_finds_the_dispersive_eigenvalues() {
    let (p, all) = drude(40, 1.0, 0.2);
    // the real parts of the modes sit near ±sqrt(μ + ω_p²) with small imaginary parts
    let region = Ellipse { center: c(1.3, -0.06), a: 0.2, b: 0.15 };
    let want = inside(&all, &region);
    assert!(want.len() >= 6, "{}", want.len());
    let r = beyn_contour(&p, &BeynOptions { nodes: 64, ..BeynOptions::new(region) }).unwrap();
    assert_same_sets(&r.values, &want, 1e-6);
    assert!(r.residuals.iter().all(|&x| x < 1e-7), "{:?}", r.residuals);
}

#[test]
fn undersized_probing_is_flagged_as_saturated() {
    let (p, all) = damped(60, 0.3, 0.02);
    let region = Ellipse { center: c(-0.15, 0.4), a: 0.1, b: 0.165 };
    let n_in = inside(&all, &region).len();
    let opts = BeynOptions { block: 2, moments: 2, nodes: 48, ..BeynOptions::new(region) };
    assert!(2 * 2 < n_in, "{n_in}");
    let r = beyn_contour(&p, &opts).unwrap();
    assert!(r.saturated, "rank {} of {}", r.rank, 4);
    assert!(r.values.len() < n_in);
}

#[test]
fn residual_inverse_iteration_converges_to_a_nearby_root() {
    for (p, all, sigma) in [
        {
            let (p, a) = drude(40, 1.0, 0.2);
            (p, a, c(1.35, -0.05))
        },
        {
            let (p, a) = damped(60, 0.3, 0.02);
            (p, a, c(-0.14, 0.52))
        },
    ] {
        let r = residual_inverse_iteration(&p, sigma, None, &RiiOptions::default()).unwrap();
        assert!(r.converged, "residual {} after {} iterations", r.residual, r.iterations);
        let d = all.iter().map(|z| (*z - r.value).norm()).fold(f64::INFINITY, f64::min);
        assert!(d < 1e-9, "{} is {d:e} from the nearest true root", r.value);
        assert!(r.residual < 1e-11);
        assert!((norm(&r.vector) - 1.0).abs() < 1e-12);
    }
}

#[test]
fn rii_polishes_a_beyn_pair() {
    let (p, _all) = drude(40, 1.0, 0.2);
    let region = Ellipse { center: c(1.3, -0.06), a: 0.2, b: 0.15 };
    let b = beyn_contour(&p, &BeynOptions { nodes: 64, ..BeynOptions::new(region) }).unwrap();
    assert!(!b.values.is_empty());
    let i = 0;
    // shift slightly off the eigenvalue (a shift on it makes T(σ) singular)
    let r = residual_inverse_iteration(
        &p,
        b.values[i] + c(2e-3, 1e-3),
        Some(b.vectors.col(i)),
        &RiiOptions { lambda0: Some(b.values[i]), ..Default::default() },
    )
    .unwrap();
    assert!(r.residual <= b.residuals[i] && r.residual < 1e-10, "{} -> {}", b.residuals[i], r.residual);
}

#[test]
fn invalid_inputs() {
    assert!(NonlinearEigenProblem::new(vec![]).is_err());
    let (p, _) = damped(20, 0.3, 0.02);
    let region = Ellipse::circle(c(0.0, 0.0), 1.0);
    assert!(beyn_contour(&p, &BeynOptions { block: 0, ..BeynOptions::new(region) }).is_err());
    assert!(beyn_contour(&p, &BeynOptions { block: 8, moments: 8, ..BeynOptions::new(region) }).is_err());
    assert!(residual_inverse_iteration(&p, c(0.0, 0.5), Some(&[C64::ONE; 3]), &RiiOptions::default()).is_err());
}
