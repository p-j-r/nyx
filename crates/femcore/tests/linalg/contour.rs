// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::contour::{ContourOptions, Ellipse, contour_eigs};
use femcore::linalg::csr::{laplacian_1d, laplacian_2d, mass_1d};
use femcore::linalg::{C64, Csr};

fn fe_spectrum(n: usize) -> Vec<f64> {
    (1..=n)
        .map(|k| {
            let c = (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos();
            6.0 * (1.0 - c) / (2.0 + c)
        })
        .collect()
}

fn sorted_re(v: &[C64]) -> Vec<C64> {
    let mut v = v.to_vec();
    v.sort_by(|a, b| (a.re, a.im).partial_cmp(&(b.re, b.im)).unwrap());
    v
}

#[test]
fn real_generalised_problem_all_eigenvalues_in_a_window() {
    let n = 120;
    let (a, m) = (laplacian_1d(n), mass_1d(n, 1.0));
    let ev = fe_spectrum(n);
    // modes 20..=29 (0-based 19..=28); the boundary tips sit in the gaps next to them
    let (lo, hi) = (0.5 * (ev[18] + ev[19]), 0.5 * (ev[28] + ev[29]));
    let opts =
        ContourOptions { nodes: 48, block: 12, moments: 6, ..ContourOptions::new(Ellipse::real_interval(lo, hi, 0.6)) };
    let r = contour_eigs(&a, Some(&m), &opts).unwrap();
    assert_eq!(r.values.len(), 10, "values {:?} (rank {}, rejected {})", r.values, r.rank, r.rejected);
    for (i, lam) in r.values.iter().enumerate() {
        assert!((lam.re - ev[19 + i]).abs() < 1e-9 && lam.im.abs() < 1e-9, "{i}: {lam} vs {}", ev[19 + i]);
    }
    assert!(r.residuals.iter().all(|&x| x < 1e-9), "{:?}", r.residuals);
    assert!(
        (r.estimated_count - 10.0).abs() < 4.0 * r.count_stderr.max(0.3),
        "count estimate {} ± {}",
        r.estimated_count,
        r.count_stderr
    );
    assert!(r.rank >= 10, "rank {} cannot be below the number of eigenvalues inside", r.rank);
    assert!(!r.saturated, "rank {}", r.rank);
    // Eigenvectors are unit and satisfy the equation.
    for i in 0..10 {
        let x = r.vectors.col(i);
        assert!((norm(x) - 1.0).abs() < 1e-12);
        assert!(true_residual(&a.to_complex(), Some(&m.to_complex()), r.values[i], x) < 1e-8);
    }
}

#[test]
fn complex_symmetric_problem_with_eigenvalues_off_the_real_axis() {
    // A = L + iγ I: eigenvalues λ_k(L) + iγ on a horizontal line; a circle around part of it.
    let (n, gamma) = (300usize, 0.35);
    let a = Csr::axpby(C64::ONE, &laplacian_1d(n).to_complex(), C64::new(0.0, gamma), &Csr::identity(n));
    let ev = laplacian_1d_eigenvalues(n);
    let (i0, i1) = (100usize, 107usize); // eight eigenvalues
    let (lo, hi) = (0.5 * (ev[i0 - 1] + ev[i0]), 0.5 * (ev[i1] + ev[i1 + 1]));
    let region = Ellipse { center: C64::new(0.5 * (lo + hi), gamma), a: 0.5 * (hi - lo), b: 0.5 * (hi - lo) };
    let opts = ContourOptions { nodes: 40, block: 10, moments: 4, ..ContourOptions::new(region) };
    let r = contour_eigs(&a, None, &opts).unwrap();
    assert_eq!(r.values.len(), 8, "{:?} rank {} rejected {}", r.values, r.rank, r.rejected);
    for (k, lam) in sorted_re(&r.values).iter().enumerate() {
        assert!((*lam - C64::new(ev[i0 + k], gamma)).norm() < 1e-9, "{k}: {lam}");
    }
    assert!(
        (r.estimated_count - 8.0).abs() < 4.0 * r.count_stderr.max(0.3),
        "{} ± {}",
        r.estimated_count,
        r.count_stderr
    );
    assert!(r.rank >= 8);
}

#[test]
fn complex_pencil_matches_a_dense_reference() {
    // Absorbing-layer pencil A x = λ B x with complex diagonal B on a small grid.
    let (nx, ny) = (20usize, 16usize);
    let n = nx * ny;
    let a = laplacian_2d(nx, ny).to_complex();
    let t: Vec<_> = (0..n)
        .map(|p| {
            let (i, j) = (p % nx, p / nx);
            let d = i.min(nx - 1 - i).min(j.min(ny - 1 - j)) as f64;
            let sg = if d < 3.0 { 2.0 * ((3.0 - d) / 3.0).powi(2) } else { 0.0 };
            (p, p, C64::new(1.0, -sg))
        })
        .collect();
    let b = Csr::from_triplets(n, n, &t).unwrap();
    let da = faer::Mat::<faer::c64>::from_fn(n, n, |i, j| a.to_dense()[i * n + j]);
    let db = faer::Mat::<faer::c64>::from_fn(n, n, |i, j| b.to_dense()[i * n + j]);
    let g = da.generalized_eigen(&db).unwrap();
    let mut ref_vals: Vec<C64> = (0..n).map(|i| g.S_a().column_vector()[i] / g.S_b().column_vector()[i]).collect();
    ref_vals.sort_by(|x, y| (x.re, x.im).partial_cmp(&(y.re, y.im)).unwrap());
    // Circle around the eigenvalue nearest to 1.2, with the radius at the largest gap in the
    // sorted distances among the 3rd..10th neighbours.
    let c0 = *ref_vals
        .iter()
        .min_by(|x, y| (**x - C64::new(1.2, 0.0)).norm().total_cmp(&(**y - C64::new(1.2, 0.0)).norm()))
        .unwrap();
    let mut d: Vec<f64> = ref_vals.iter().map(|z| (*z - c0).norm()).collect();
    d.sort_by(f64::total_cmp);
    let k = (3..10).max_by(|&p, &q| (d[p + 1] / d[p]).total_cmp(&(d[q + 1] / d[q]))).unwrap();
    let radius = 0.5 * (d[k] + d[k + 1]);
    let inside: Vec<C64> = ref_vals.iter().copied().filter(|z| (*z - c0).norm() <= radius).collect();
    let opts = ContourOptions { nodes: 48, block: 10, moments: 5, ..ContourOptions::new(Ellipse::circle(c0, radius)) };
    let r = contour_eigs(&a, Some(&b), &opts).unwrap();
    assert_eq!(r.values.len(), inside.len(), "found {:?}\nwant  {:?}", r.values, inside);
    for (g, w) in sorted_re(&r.values).iter().zip(&sorted_re(&inside)) {
        assert!((*g - *w).norm() < 1e-8 * w.norm().max(1.0), "{g} vs {w}");
    }
    assert!(r.residuals.iter().all(|&x| x < 1e-8));
}

#[test]
fn saturation_flag_and_count_estimate_warn_when_the_subspace_is_too_small() {
    let n = 120;
    let (a, m) = (laplacian_1d(n), mass_1d(n, 1.0));
    let ev = fe_spectrum(n);
    let (lo, hi) = (0.5 * (ev[18] + ev[19]), 0.5 * (ev[28] + ev[29])); // 10 eigenvalues
    let opts =
        ContourOptions { nodes: 48, block: 2, moments: 3, ..ContourOptions::new(Ellipse::real_interval(lo, hi, 0.6)) };
    let r = contour_eigs(&a, Some(&m), &opts).unwrap();
    assert!(r.saturated, "rank {}", r.rank);
    assert!(r.values.len() < 10);
    // The stochastic estimate still sees all of them (up to its variance), flagging the shortfall.
    assert!(
        r.estimated_count + 3.0 * r.count_stderr > r.values.len() as f64 + 1.0,
        "estimate {} ± {} found {}",
        r.estimated_count,
        r.count_stderr,
        r.values.len()
    );
}

#[test]
fn empty_region_returns_nothing() {
    let n = 120;
    let (a, m) = (laplacian_1d(n), mass_1d(n, 1.0));
    let ev = fe_spectrum(n);
    let (lo, hi) = (ev[40] + 0.2 * (ev[41] - ev[40]), ev[40] + 0.8 * (ev[41] - ev[40])); // inside one gap
    let opts = ContourOptions { nodes: 32, ..ContourOptions::new(Ellipse::real_interval(lo, hi, 1.0)) };
    let r = contour_eigs(&a, Some(&m), &opts).unwrap();
    assert!(r.values.is_empty(), "{:?}", r.values);
    assert!(
        r.estimated_count.abs() < 4.0 * r.count_stderr.max(0.1) + 0.2,
        "{} ± {}",
        r.estimated_count,
        r.count_stderr
    );
}

#[test]
fn results_do_not_depend_on_the_number_of_workers() {
    let n = 200;
    let (a, m) = (laplacian_1d(n), mass_1d(n, 1.0));
    let ev = fe_spectrum(n);
    let (lo, hi) = (0.5 * (ev[28] + ev[29]), 0.5 * (ev[35] + ev[36]));
    let mk = |workers| ContourOptions {
        nodes: 24,
        block: 10,
        moments: 4,
        workers,
        ..ContourOptions::new(Ellipse::real_interval(lo, hi, 0.7))
    };
    let one = contour_eigs(&a, Some(&m), &mk(1)).unwrap();
    let three = contour_eigs(&a, Some(&m), &mk(3)).unwrap();
    assert_eq!(one.values.len(), 7);
    assert_eq!(one.values, three.values, "eigenvalues differ between 1 and 3 workers");
    assert_eq!(one.rank, three.rank);
}

#[test]
fn invalid_arguments_are_rejected() {
    let a = laplacian_1d(30);
    let e = Ellipse::circle(C64::ONE, 0.1);
    assert!(contour_eigs(&a, None, &ContourOptions { nodes: 2, ..ContourOptions::new(e) }).is_err());
    assert!(contour_eigs(&a, None, &ContourOptions { block: 30, moments: 4, ..ContourOptions::new(e) }).is_err());
    assert!(contour_eigs(&a, None, &ContourOptions::new(Ellipse { center: C64::ONE, a: 0.0, b: 1.0 })).is_err());
}
