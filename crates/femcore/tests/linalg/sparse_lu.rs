// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, laplacian_2d, laplacian_3d, mass_1d};
use femcore::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};
use femcore::linalg::rng::random_vector;
use femcore::linalg::sparse_lu::{LuOptions, LuOrdering, ShiftedSparseLu, SparseLu};
use femcore::linalg::{C64, Csr, Scalar};

fn rel_residual<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::ZERO; b.len()];
    a.spmv(x, &mut ax);
    ax.iter().zip(b).map(|(p, q)| (*q - *p).modulus_sq()).sum::<f64>().sqrt() / norm(b)
}

fn helmholtz(nx: usize, k2: f64, sigma: f64) -> Csr<C64> {
    let l = laplacian_2d(nx, nx).to_complex();
    let t: Vec<_> = (0..nx * nx).map(|p| (p, p, C64::new(-k2, sigma * (1.0 + (p % nx) as f64 / nx as f64)))).collect();
    Csr::axpby(C64::ONE, &l, C64::ONE, &Csr::from_triplets(nx * nx, nx * nx, &t).unwrap())
}

fn all_orderings(n: usize) -> Vec<(&'static str, LuOrdering)> {
    vec![
        ("natural", LuOrdering::Natural),
        ("amd", LuOrdering::Amd),
        ("nd", LuOrdering::NestedDissection { leaf: 16 }),
        ("reverse", LuOrdering::Custom((0..n).rev().collect())),
    ]
}

#[test]
fn every_ordering_solves_real_and_complex_systems() {
    let a = laplacian_2d(24, 24);
    let b = random_vector::<f64>(a.nrows(), 1);
    for (name, ordering) in all_orderings(a.nrows()) {
        let lu = SparseLu::factor(&a, &LuOptions { ordering, pivot_tol: 0.1 }).unwrap();
        let mut x = b.clone();
        lu.solve_in_place(&mut x);
        assert!(rel_residual(&a, &x, &b) < 1e-12, "{name}: {}", rel_residual(&a, &x, &b));
    }
    let h = helmholtz(24, 1.3, 0.2);
    let bc = random_vector::<C64>(h.nrows(), 2);
    for (name, ordering) in all_orderings(h.nrows()) {
        let lu = SparseLu::factor(&h, &LuOptions { ordering, pivot_tol: 0.1 }).unwrap();
        let mut x = bc.clone();
        lu.solve_in_place(&mut x);
        assert!(rel_residual(&h, &x, &bc) < 1e-11, "{name}: {}", rel_residual(&h, &x, &bc));
    }
}

#[test]
fn pivoting_handles_zero_diagonals_and_nonsymmetric_patterns() {
    // a randomly row-permuted well-conditioned matrix has zeros on the diagonal
    let n = 120;
    let base = laplacian_1d(n);
    let perm: Vec<usize> = (0..n).map(|i| (i * 37 + 11) % n).collect();
    let mut t = Vec::new();
    for i in 0..n {
        let (cols, vals) = base.row(i);
        for (&c, &v) in cols.iter().zip(vals) {
            t.push((perm[i], c as usize, v));
        }
    }
    let a = Csr::from_triplets(n, n, &t).unwrap();
    assert!(a.diagonal().contains(&0.0));
    let b = random_vector::<f64>(n, 3);
    for tol in [1.0, 0.1, 0.001] {
        for ordering in [LuOrdering::Natural, LuOrdering::Amd] {
            let lu = SparseLu::factor(&a, &LuOptions { ordering, pivot_tol: tol }).unwrap();
            assert!(lu.stats().off_diagonal_pivots > 0);
            let mut x = b.clone();
            lu.solve_in_place(&mut x);
            assert!(rel_residual(&a, &x, &b) < 1e-10, "tol {tol}: {}", rel_residual(&a, &x, &b));
        }
    }
}

#[test]
fn singular_matrices_and_bad_orderings_are_reported() {
    let a = Csr::from_triplets(3, 3, &[(0, 0, 1.0), (1, 1, 1.0), (2, 0, 1.0)]).unwrap(); // column 2 empty
    assert_eq!(SparseLu::factor(&a, &LuOptions::default()).unwrap_err(), SolverError::Singular);
    let ok = laplacian_1d(4);
    let bad = LuOptions { ordering: LuOrdering::Custom(vec![0, 1, 1, 3]), pivot_tol: 0.1 };
    assert!(matches!(SparseLu::factor(&ok, &bad), Err(SolverError::InvalidArgument(_))));
}

#[test]
fn ordering_reduces_fill_and_flops_on_a_3d_grid() {
    let a = laplacian_3d(14);
    let nat = SparseLu::factor(&a, &LuOptions { ordering: LuOrdering::Natural, pivot_tol: 0.1 }).unwrap();
    let amd = SparseLu::factor(&a, &LuOptions { ordering: LuOrdering::Amd, pivot_tol: 0.1 }).unwrap();
    let nd = SparseLu::factor(&a, &LuOptions { ordering: LuOrdering::NestedDissection { leaf: 50 }, pivot_tol: 0.1 })
        .unwrap();
    println!(
        "3-D 14³: nnz(L+U) natural {} / amd {} / nd {}; flops {:.2e} / {:.2e} / {:.2e}",
        nat.factor_nnz(),
        amd.factor_nnz(),
        nd.factor_nnz(),
        nat.stats().flops as f64,
        amd.stats().flops as f64,
        nd.stats().flops as f64
    );
    assert!(amd.factor_nnz() * 2 < nat.factor_nnz());
    assert!(nd.factor_nnz() * 2 < nat.factor_nnz());
    assert!(amd.stats().flops * 4 < nat.stats().flops);
    // for a symmetric pattern the diagonal is always accepted on an SPD matrix
    assert_eq!(nd.stats().off_diagonal_pivots, 0);
}

#[test]
fn shifted_wrapper_agrees_with_the_faer_factorisation_on_a_pencil() {
    let n = 200;
    let (a, m) = (laplacian_1d(n), mass_1d(n, 1.0 / (n as f64 + 1.0)));
    let sigma = 3.7;
    let b = random_vector::<f64>(n, 4);
    let mut reference = ShiftedFactorization::new(&a, Some(&m), Definiteness::General).unwrap();
    reference.set_shift(sigma).unwrap();
    let mut x_ref = b.clone();
    reference.solve_block_in_place(&mut x_ref, 1).unwrap();
    for ordering in [LuOrdering::Amd, LuOrdering::NestedDissection { leaf: 8 }] {
        let mut ours = ShiftedSparseLu::new(&a, Some(&m), &LuOptions { ordering, pivot_tol: 0.1 }).unwrap();
        ours.set_shift(sigma).unwrap();
        let mut x = b.clone();
        ours.solve_block_in_place(&mut x, 1).unwrap();
        let err = x.iter().zip(&x_ref).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&x_ref);
        assert!(err < 1e-8, "{err}");
        // a second shift reuses the ordering
        ours.set_shift(1.1).unwrap();
        assert_eq!(ours.stats().numeric_count, 2);
        assert_eq!(ours.stats().symbolic_count, 1);
    }
}

#[test]
fn block_solve_matches_single_solves() {
    let a = helmholtz(12, 0.8, 0.3);
    let lu = SparseLu::factor(&a, &LuOptions::default()).unwrap();
    let n = a.nrows();
    let rhs = random_vector::<C64>(3 * n, 6);
    let mut blk = rhs.clone();
    lu.solve_block_in_place(&mut blk, 3);
    for c in 0..3 {
        let mut one = rhs[c * n..(c + 1) * n].to_vec();
        lu.solve_in_place(&mut one);
        assert_eq!(&blk[c * n..(c + 1) * n], &one[..]);
    }
}
