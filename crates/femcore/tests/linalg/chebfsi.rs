// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Chebyshev-filtered subspace iteration: many lowest eigenpairs of a 2-D Laplacian against the closed form,
//! complex Hermitian operators, the effect of the degree, and input validation.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::chebfsi::{ChebOptions, chebyshev_subspace};
use femcore::linalg::csr::{laplacian_1d, laplacian_2d};
use femcore::linalg::{C64, Csr, Scalar};

#[test]
fn forty_lowest_modes_of_the_2d_laplacian() {
    let (nx, ny) = (50, 50);
    let a = laplacian_2d(nx, ny);
    let exact = laplacian_2d_eigenvalues(nx, ny);
    let r = chebyshev_subspace(&a, &ChebOptions { tol: 1e-8, max_it: 200, ..ChebOptions::new(40) }).unwrap();
    assert!(r.converged, "{} iterations, locked {:?}", r.iterations, r.locked_history.last());
    for (g, w) in r.values.iter().zip(&exact) {
        assert!((g - w).abs() < 1e-7 * w, "{g} vs {w}");
    }
    assert!(r.residuals.iter().all(|&e| e < 1e-6), "{:?}", r.residuals.iter().cloned().fold(0.0, f64::max));
    // orthonormal vectors
    for i in 0..40 {
        for j in 0..=i {
            let d: f64 = r.vectors[i].iter().zip(&r.vectors[j]).map(|(p, q)| p * q).sum();
            assert!((d - if i == j { 1.0 } else { 0.0 }).abs() < 1e-6, "({i},{j}) {d}");
        }
    }
    // locking is monotone and happens progressively
    assert!(r.locked_history.windows(2).all(|w| w[1] >= w[0]));
    eprintln!(
        "Chebyshev FSI: 40 modes of the {nx}x{ny} Laplacian in {} iterations, {} matvecs",
        r.iterations, r.matvecs
    );
}

#[test]
fn degree_trades_iterations_against_matrix_products() {
    let a = laplacian_2d(40, 40);
    let mut rows = Vec::new();
    for degree in [8usize, 20, 40] {
        let r =
            chebyshev_subspace(&a, &ChebOptions { tol: 1e-8, degree, max_it: 400, ..ChebOptions::new(20) }).unwrap();
        assert!(r.converged, "degree {degree}");
        rows.push((degree, r.iterations, r.matvecs));
    }
    eprintln!("degree / iterations / matvecs: {rows:?}");
    assert!(rows[0].1 > rows[2].1, "a higher degree should need fewer outer iterations: {rows:?}");
}

#[test]
fn complex_hermitian_operator() {
    // phase-twisted tridiagonal chain: spectrum of the 1-D Laplacian
    let m = 200;
    let mut t = Vec::new();
    for i in 0..m {
        t.push((i, i, C64::new(2.0, 0.0)));
        if i + 1 < m {
            let ph = C64::new(0.0, 0.2 * (1 + i % 4) as f64).exp();
            t.push((i, i + 1, ph.scale(-1.0)));
            t.push((i + 1, i, ph.conj().scale(-1.0)));
        }
    }
    let h = Csr::from_triplets(m, m, &t).unwrap();
    let r =
        chebyshev_subspace(&h, &ChebOptions { tol: 1e-8, degree: 40, max_it: 400, ..ChebOptions::new(10) }).unwrap();
    assert!(
        r.converged,
        "{} iterations, locked {:?}, residuals {:?}",
        r.iterations,
        r.locked_history.last(),
        r.residuals
    );
    for (g, w) in r.values.iter().zip(laplacian_1d_eigenvalues(m)) {
        assert!((g - w).abs() < 1e-7, "{g} vs {w}");
    }
}

#[test]
fn invalid_inputs() {
    let a = laplacian_1d(30);
    assert!(chebyshev_subspace(&a, &ChebOptions::new(0)).is_err());
    assert!(chebyshev_subspace(&a, &ChebOptions::new(14)).is_err());
    assert!(chebyshev_subspace(&a, &ChebOptions { degree: 1, ..ChebOptions::new(3) }).is_err());
    let _ = Scalar::modulus(1.0f64);
}
