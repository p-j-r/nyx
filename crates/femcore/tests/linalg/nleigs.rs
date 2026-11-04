#![allow(clippy::needless_range_loop)]
//! NLEIGS against the closed-form spectra of the scalar-reducible problems of `nep.rs`, and against Beyn's method.
use crate::common;

use common::*;
use femcore::linalg::contour::Ellipse;
use femcore::linalg::csr::laplacian_1d;
use femcore::linalg::dense::{Dense, solve_small};
use femcore::linalg::faer_glue::generalized_eig_c64;
use femcore::linalg::linop::LinearOperator;
use femcore::linalg::nep::{BeynOptions, NepTerm, NonlinearEigenProblem, beyn_contour};
use femcore::linalg::nleigs::{NleigsOperator, NleigsOptions, build_interpolant, nleigs, ray_samples};
use femcore::linalg::rng::random_vector;
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

/// `L − λ I + i √λ I`: radiation-type condition; scalar equation `μ − λ + i√λ = 0`, `√λ = (i + √(4μ−1))/2`,
/// so `λ = μ − ½ + i √(4μ−1)/2` for `μ > ¼` (principal branch).
fn sqrt_problem(n: usize) -> (NonlinearEigenProblem, Vec<C64>) {
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let p = NonlinearEigenProblem::new(vec![
        NepTerm { matrix: l, f: Box::new(|_| C64::ONE) },
        NepTerm { matrix: id, f: Box::new(|z| -z + C64::new(0.0, 1.0) * z.sqrt()) },
    ])
    .unwrap();
    let all = laplacian_1d_eigenvalues(n)
        .into_iter()
        .filter(|&mu| mu > 0.25)
        .map(|mu| c(mu - 0.5, 0.5 * (4.0 * mu - 1.0).sqrt()))
        .collect();
    (p, all)
}

#[test]
fn rational_interpolation_of_entire_branch_and_pole_functions() {
    let id: Csr<C64> = Csr::identity(2);
    let mk = |f: Box<dyn Fn(C64) -> C64>| NonlinearEigenProblem::new(vec![NepTerm { matrix: id.clone(), f }]).unwrap();
    let disc = Ellipse::circle(c(3.0, 0.0), 2.5);
    // entire: polynomial (all poles at infinity) interpolation converges geometrically
    let (it, err) = build_interpolant(&mk(Box::new(|z| (-z).exp())), &disc, &[], 1e-12, 80, 400).unwrap();
    assert!(err < 1e-11, "exp: {err} at degree {}", it.degree());
    // a pole outside the region: the rational basis handles it with a modest degree, a polynomial would need far more
    let pole = c(-1.0, 0.5);
    let f_pole = move |z: C64| C64::ONE / (z - pole);
    let (it, err) = build_interpolant(&mk(Box::new(f_pole)), &disc, &[pole], 1e-12, 80, 400).unwrap();
    assert!(err < 1e-10 && it.degree() <= 24, "pole: {err} at degree {}", it.degree());
    let (itp, errp) = build_interpolant(&mk(Box::new(f_pole)), &disc, &[], 1e-12, 24, 400).unwrap();
    assert!(errp > 100.0 * err, "polynomial {errp} vs rational {err} (degrees {} / {})", itp.degree(), it.degree());
    // branch point at 0 with the cut along the negative axis
    let cut = ray_samples(C64::ZERO, c(-1.0, 0.0), 1e-3, 1e4, 60);
    let (it, err) = build_interpolant(&mk(Box::new(|z| z.sqrt())), &disc, &cut, 1e-10, 80, 400).unwrap();
    assert!(err < 1e-8, "sqrt: {err} at degree {}", it.degree());
}

#[test]
fn the_dense_linearisation_has_the_nonlinear_eigenvalues_and_the_operator_matches_it() {
    let (p, all) = sqrt_problem(3);
    let region = Ellipse::circle(c(1.6, 1.2), 0.9);
    let cut = ray_samples(C64::ZERO, c(-1.0, 0.0), 1e-3, 1e4, 60);
    let (it, _) = build_interpolant(&p, &region, &cut, 1e-12, 60, 300).unwrap();
    let s = region.center + c(0.3, 0.2); // a shift away from the centre exercises the relative arithmetic
    let op = NleigsOperator::new(&p, it, s).unwrap();
    let (a, b) = op.dense_pencil();
    // eigenvalues of the pencil (relative to the region's centre) inside the region are the true ones
    let (vals, _) = generalized_eig_c64(&b, &a).unwrap();
    let want: Vec<C64> = all.iter().copied().filter(|z| region.contains(*z)).collect();
    assert!(!want.is_empty(), "test region contains no eigenvalue");
    for w in &want {
        let d = vals.iter().map(|v| (*v + region.center - *w).norm()).fold(f64::INFINITY, f64::min);
        assert!(d < 1e-8, "eigenvalue {w} missing from the pencil (closest {d:e})");
    }
    // the matrix-free operator equals (B − sA)⁻¹ A v
    let n = op.dim();
    let v = random_vector::<C64>(n, 3);
    let mut y = vec![C64::ZERO; n];
    op.apply(&v, &mut y);
    let av: Vec<C64> = (0..n).map(|i| (0..n).fold(C64::ZERO, |sum, j| sum + a[(i, j)] * v[j])).collect();
    let sp = s - region.center; // the pencil is expressed relative to the centre
    let m = Dense::from_fn(n, n, |i, j| b[(i, j)] - sp * a[(i, j)]);
    let want_y = solve_small(&m, &av).unwrap();
    let err = y.iter().zip(&want_y).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm(&want_y);
    assert!(err < 1e-9, "operator vs dense: {err}");
}

#[test]
fn nleigs_finds_all_quadratic_eigenvalues_in_a_region() {
    let (p, all) = damped(60, 0.3, 0.02);
    let region = Ellipse { center: c(-0.15, 0.4), a: 0.12, b: 0.18 };
    let want = inside(&all, &region);
    assert!(want.len() >= 6, "{}", want.len());
    let r = nleigs(&p, &NleigsOptions { nev: want.len() + 6, ..NleigsOptions::new(region, 0) }).unwrap();
    assert_same_sets(&r.values, &want, 1e-8);
    assert!(r.residuals.iter().all(|&x| x < 1e-10), "{:?}", r.residuals);
    println!("damped: degree {}, {} matvecs, {} restarts", r.degree, r.matvecs, r.restarts);
}

#[test]
fn nleigs_handles_the_pole_of_a_dispersive_medium() {
    let (p, all) = drude(40, 1.0, 0.2);
    let region = Ellipse { center: c(1.3, -0.06), a: 0.2, b: 0.15 };
    let want = inside(&all, &region);
    assert!(want.len() >= 6);
    // the only singularity of the f_t is the pole at −iγ
    let opts =
        NleigsOptions { singularities: vec![c(0.0, -0.2)], nev: want.len() + 6, ..NleigsOptions::new(region, 0) };
    let r = nleigs(&p, &opts).unwrap();
    assert_same_sets(&r.values, &want, 1e-8);
    assert!(r.residuals.iter().all(|&x| x < 1e-10), "{:?}", r.residuals);
    assert!(r.interpolation_error < 1e-9);
    println!(
        "drude: degree {}, interpolation error {:.1e}, max residual {:.1e}",
        r.degree,
        r.interpolation_error,
        r.residuals.iter().cloned().fold(0.0, f64::max)
    );
}

#[test]
fn nleigs_handles_a_square_root_branch_cut() {
    let (p, all) = sqrt_problem(50);
    let region = Ellipse { center: c(2.5, 1.0), a: 2.0, b: 1.2 };
    let want = inside(&all, &region);
    assert!(want.len() >= 5, "{}", want.len());
    let opts = NleigsOptions {
        singularities: ray_samples(C64::ZERO, c(-1.0, 0.0), 1e-3, 1e4, 80),
        nev: want.len() + 8,
        ..NleigsOptions::new(region, 0)
    };
    let r = nleigs(&p, &opts).unwrap();
    println!(
        "sqrt: degree {}, interpolation error {:.1e}, {} found of {}",
        r.degree,
        r.interpolation_error,
        r.values.len(),
        want.len()
    );
    assert_same_sets(&r.values, &want, 1e-7);
    assert!(r.residuals.iter().all(|&x| x < 1e-9), "{:?}", r.residuals);
}

#[test]
fn nleigs_and_beyn_agree_on_a_delay_problem() {
    // T(λ) = L − λ I + e^{−λ} I : entire, no closed form
    let n = 30;
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let p = NonlinearEigenProblem::new(vec![
        NepTerm { matrix: l, f: Box::new(|_| C64::ONE) },
        NepTerm { matrix: id, f: Box::new(|z| -z + (-z).exp()) },
    ])
    .unwrap();
    let region = Ellipse { center: c(1.0, 0.0), a: 1.5, b: 1.0 };
    let beyn = beyn_contour(&p, &BeynOptions { nodes: 64, ..BeynOptions::new(region) }).unwrap();
    assert!(beyn.values.len() >= 3, "{:?}", beyn.values);
    let r = nleigs(&p, &NleigsOptions { nev: beyn.values.len() + 6, ..NleigsOptions::new(region, 0) }).unwrap();
    assert_same_sets(&r.values, &beyn.values, 1e-7);
}
