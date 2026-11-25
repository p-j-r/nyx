// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Differential verification for `femcore` linear algebra kernels across default (SIMD/`faer`)
//! and `--features portable-kernels` configurations.
//!
//! Runs identically under both feature flags (`cargo test --release -p femcore --test portable_kernels`
//! and `cargo test --release -p femcore --features portable-kernels --test portable_kernels`) and
//! checks every solve, operator product, and eigenproblem against exact/reference residuals `< 1e-13`.

use femcore::linalg::basis::Basis;
use femcore::linalg::csr::{laplacian_1d, laplacian_2d, Csr};
use femcore::linalg::direct::FactorBackend;
use femcore::linalg::driver::{
    dense_eigh, dense_gen_eigh, eigs_near, eigs_real_general_near, DirectSolver, EigenOptions,
};
use femcore::linalg::faer_glue::generalized_eig_c64;
use femcore::linalg::dense::Dense;
use femcore::linalg::problem::Symmetry;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::rng::SplitMix64;
use femcore::linalg::scalar::{Scalar, C64};
use femcore::linalg::vecops;

fn rand_vec<T: Scalar>(n: usize, seed: u64) -> Vec<T> {
    let mut r = SplitMix64::new(seed);
    (0..n).map(|_| T::from_parts(r.next_signed(), r.next_signed())).collect()
}

fn vec_norm2<T: Scalar>(v: &[T]) -> f64 {
    v.iter().map(|x| x.modulus_sq()).sum::<f64>().sqrt()
}

fn rel_res<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::ZERO; a.nrows()];
    a.spmv(x, &mut ax);
    let diff = ax.iter().zip(b).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt();
    diff / vec_norm2(b).max(1e-300)
}

#[test]
fn direct_solver_lu_and_lblt_achieve_sub_1e13_residual() {
    // 1. Real unsymmetric sparse LU (`LuColamd`).
    let n = 48;
    let mut unsym_f64 = Vec::new();
    for i in 0..n {
        unsym_f64.push((i, i, 5.0 + 0.05 * i as f64));
        if i + 1 < n {
            unsym_f64.push((i, i + 1, -1.2));
            unsym_f64.push((i + 1, i, -0.4));
        }
        if i + 3 < n {
            unsym_f64.push((i, i + 3, 0.15));
        }
    }
    let a_lu_f64 = Csr::from_triplets(n, n, &unsym_f64).unwrap();
    let solver_lu = DirectSolver::new(&a_lu_f64, false).unwrap();
    assert_eq!(solver_lu.stats().backend, FactorBackend::LuColamd);
    let b0 = rand_vec::<f64>(n, 101);
    let b1 = rand_vec::<f64>(n, 102);
    let x0 = solver_lu.solve(&b0);
    assert!(rel_res(&a_lu_f64, &x0, &b0) < 1e-13, "real LU single RHS residual");
    let xs = solver_lu.solve_many(&[b0.clone(), b1.clone()]);
    assert!(rel_res(&a_lu_f64, &xs[0], &b0) < 1e-13);
    assert!(rel_res(&a_lu_f64, &xs[1], &b1) < 1e-13);

    // 2. Complex unsymmetric sparse LU (`LuColamd`), plus `update_matrix` symbolic reuse.
    let mut unsym_c64 = Vec::new();
    for i in 0..n {
        unsym_c64.push((i, i, C64::new(4.5, 0.25)));
        if i + 1 < n {
            unsym_c64.push((i, i + 1, C64::new(-1.0, 0.15)));
            unsym_c64.push((i + 1, i, C64::new(-0.35, -0.25)));
        }
    }
    let a_lu_c64 = Csr::from_triplets(n, n, &unsym_c64).unwrap();
    let mut solver_lu_c = DirectSolver::new(&a_lu_c64, false).unwrap();
    assert_eq!(solver_lu_c.stats().backend, FactorBackend::LuColamd);
    let bc0 = rand_vec::<C64>(n, 201);
    let xc0 = solver_lu_c.solve(&bc0);
    assert!(rel_res(&a_lu_c64, &xc0, &bc0) < 1e-13, "complex LU residual");

    // Symbolic reuse on updated numeric values with the exact same sparsity pattern.
    let mut unsym_c64_2 = unsym_c64;
    for (r, c, v) in &mut unsym_c64_2 {
        if r == c {
            *v = C64::new(6.0, -0.4);
        }
    }
    let a_lu_c64_2 = Csr::from_triplets(n, n, &unsym_c64_2).unwrap();
    solver_lu_c.update_matrix(&a_lu_c64_2).unwrap();
    assert_eq!(solver_lu_c.stats().symbolic_count, 1);
    let xc1 = solver_lu_c.solve(&bc0);
    assert!(rel_res(&a_lu_c64_2, &xc1, &bc0) < 1e-13, "updated complex LU residual");

    // 3. Real symmetric indefinite LBL^T (`LbltAmd` / `LbltNd`).
    let a_indef = Csr::axpby(1.0, &laplacian_1d(n), -0.35, &Csr::identity(n));
    let solver_lblt_r = DirectSolver::new(&a_indef, false).unwrap();
    assert!(
        matches!(solver_lblt_r.stats().backend, FactorBackend::LbltAmd | FactorBackend::LbltNd),
        "expected LBLT backend for real symmetric indefinite, got {:?}",
        solver_lblt_r.stats().backend
    );
    let xr = solver_lblt_r.solve(&b0);
    assert!(rel_res(&a_indef, &xr, &b0) < 1e-13, "real LBLT residual");

    // 4. Complex-symmetric (non-Hermitian) LBL^T (`LbltAmd` / `LbltNd`).
    let mut sym_c64 = Vec::new();
    for i in 0..n {
        sym_c64.push((i, i, C64::new(2.0, -0.15) - C64::new(0.45, 0.02)));
        if i + 1 < n {
            let off = C64::new(-1.0, 0.08);
            sym_c64.push((i, i + 1, off));
            sym_c64.push((i + 1, i, off));
        }
    }
    let a_lblt_c = Csr::from_triplets(n, n, &sym_c64).unwrap();
    let solver_lblt_c = DirectSolver::new(&a_lblt_c, false).unwrap();
    assert!(
        matches!(solver_lblt_c.stats().backend, FactorBackend::LbltAmd | FactorBackend::LbltNd),
        "expected LBLT backend for complex-symmetric, got {:?}",
        solver_lblt_c.stats().backend
    );
    let x_lblt_c = solver_lblt_c.solve_many(&[bc0.clone(), rand_vec::<C64>(n, 202)]);
    assert!(rel_res(&a_lblt_c, &x_lblt_c[0], &bc0) < 1e-13, "complex LBLT residual");

    // 5. Real SPD Cholesky (`spd = true`).
    let a_spd = laplacian_2d(8, 7);
    let solver_chol = DirectSolver::new(&a_spd, true).unwrap();
    let b_spd = rand_vec::<f64>(a_spd.nrows(), 301);
    let x_chol = solver_chol.solve(&b_spd);
    assert!(rel_res(&a_spd, &x_chol, &b_spd) < 1e-13, "SPD Cholesky residual");
}

#[test]
fn sparse_ops_spmv_and_axpby_match_exact_reference() {
    let a = laplacian_2d(7, 6).to_complex();
    let b = Csr::from_diag(&rand_vec::<C64>(a.nrows(), 11));
    let alpha = C64::new(1.25, -0.5);
    let beta = C64::new(-0.75, 0.3);
    let c = Csr::axpby(alpha, &a, beta, &b);

    let x = rand_vec::<C64>(a.ncols(), 12);
    let mut y_ax = vec![C64::ZERO; a.nrows()];
    let mut y_bx = vec![C64::ZERO; a.nrows()];
    let mut y_cx = vec![C64::ZERO; a.nrows()];
    let mut y_cx_dot = vec![C64::ZERO; a.nrows()];

    a.spmv(&x, &mut y_ax);
    b.spmv(&x, &mut y_bx);
    c.spmv(&x, &mut y_cx);
    let fused_dot = c.spmv_dotc(&x, &mut y_cx_dot);
    assert_eq!(y_cx, y_cx_dot);
    assert_eq!(fused_dot, vecops::dotc(&x, &y_cx));

    let ref_y: Vec<C64> = y_ax.iter().zip(&y_bx).map(|(p, q)| alpha * *p + beta * *q).collect();
    let err = y_cx.iter().zip(&ref_y).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt();
    let rel = err / vec_norm2(&ref_y).max(1e-300);
    assert!(rel < 1e-13, "spmv/axpby relative discrepancy {rel:e}");

    let ab = a.matmul(&b);
    let mut y_abx = vec![C64::ZERO; a.nrows()];
    ab.spmv(&x, &mut y_abx);
    let mut y_a_bx = vec![C64::ZERO; a.nrows()];
    a.spmv(&y_bx, &mut y_a_bx);
    let mm_err = y_abx.iter().zip(&y_a_bx).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt();
    assert!(mm_err / vec_norm2(&y_a_bx).max(1e-300) < 1e-13);
}

#[test]
fn dense_eigensolvers_satisfy_residual_under_1e13() {
    let n = 16;
    let mut k = vec![0.0f64; n * n];
    let mut m = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let v = ((i * 19 + j * 37 + 5) % 89) as f64 / 89.0 - 0.5;
            let kv = if i == j { 8.0 + v.abs() } else { 0.15 * v };
            let mv = if i == j { 2.5 + 0.1 * v.abs() } else { 0.03 * v };
            k[i * n + j] = kv;
            k[j * n + i] = kv;
            m[i * n + j] = mv;
            m[j * n + i] = mv;
        }
    }
    let (vals, vecs) = dense_eigh(&k, n).unwrap();
    for (j, &lam) in vals.iter().enumerate() {
        let v = &vecs[j];
        let mut r2 = 0.0;
        for row in 0..n {
            let kv: f64 = (0..n).map(|col| k[row * n + col] * v[col]).sum();
            let d = kv - lam * v[row];
            r2 += d * d;
        }
        let rel = r2.sqrt() / ((lam.abs() + 1.0) * vec_norm2(v));
        assert!(rel < 1e-13, "dense_eigh mode {j} residual {rel:e}");
    }

    let (gvals, gvecs) = dense_gen_eigh(&k, &m, n).unwrap();
    for (j, &lam) in gvals.iter().enumerate() {
        let v = &gvecs[j];
        let mut r2 = 0.0;
        let mut mv2 = 0.0;
        for row in 0..n {
            let kv: f64 = (0..n).map(|col| k[row * n + col] * v[col]).sum();
            let mv: f64 = (0..n).map(|col| m[row * n + col] * v[col]).sum();
            let d = kv - lam * mv;
            r2 += d * d;
            mv2 += mv * mv;
        }
        let rel = r2.sqrt() / ((lam.abs() + 1.0) * mv2.sqrt().max(1e-300));
        assert!(rel < 1e-13, "dense_gen_eigh mode {j} residual {rel:e}");
    }

    let ac = Dense::from_fn(6, 6, |i, j| {
        C64::new((i + 2 * j) as f64 + if i == j { 4.0 } else { 0.0 }, 0.1 * (i as f64 - j as f64))
    });
    let bc = Dense::from_fn(6, 6, |i, j| {
        if i == j { C64::new(2.0 + i as f64 * 0.2, 0.15) } else { C64::new(0.05, -0.02) }
    });
    let (lam_c, vec_c) = generalized_eig_c64(&ac, &bc).unwrap();
    for col in 0..6 {
        let mut r2 = 0.0;
        let mut bx2 = 0.0;
        for i in 0..6 {
            let ax = (0..6).fold(C64::ZERO, |s, j| s + ac[(i, j)] * vec_c[(j, col)]);
            let bx = (0..6).fold(C64::ZERO, |s, j| s + bc[(i, j)] * vec_c[(j, col)]);
            r2 += (ax - lam_c[col] * bx).norm_sqr();
            bx2 += bx.norm_sqr();
        }
        let rel = r2.sqrt() / ((lam_c[col].norm() + 1.0) * bx2.sqrt().max(1e-300));
        assert!(rel < 1e-13, "generalized_eig_c64 mode {col} residual {rel:e}");
    }
}

fn check_tile_mult_and_basis<T: Scalar>() {
    for &(len, nrows, ncols, ldq) in &[(19usize, 7usize, 5usize, 9usize), (137, 15, 9, 17)] {
        let src = rand_vec::<T>(len * nrows, 31);
        let q = rand_vec::<T>(ldq * ncols, 32);
        let mut got = vec![T::ZERO; len * ncols];
        T::k_tile_mult(&src, len, nrows, &q, ldq, ncols, &mut got);

        let mut want = vec![T::ZERO; len * ncols];
        for c in 0..ncols {
            for i in 0..len {
                want[c * len + i] =
                    (0..nrows).fold(T::ZERO, |acc, j| acc + src[j * len + i] * q[c * ldq + j]);
            }
        }
        let err = got.iter().zip(&want).map(|(g, w)| (*g - *w).modulus_sq()).sum::<f64>().sqrt();
        let rel = err / vec_norm2(&want).max(1e-300);
        assert!(rel < 1e-13, "k_tile_mult ({len}x{nrows} * {nrows}x{ncols}) rel err {rel:e}");
    }

    // Exercise `Basis::mult_in_place` across multiple row tiles (`TILE_ROWS = 512`).
    let (n, nrows_q, ncols_q, ldq) = (1200usize, 7usize, 5usize, 8usize);
    let mut basis = Basis::<T>::new(n, 10);
    for j in 0..10 {
        basis.col_mut(j).copy_from_slice(&rand_vec::<T>(n, 50 + j as u64));
    }
    let orig = basis.clone();
    let q = rand_vec::<T>(ldq * ncols_q, 99);
    basis.mult_in_place(1, nrows_q, ncols_q, &q, ldq);
    for c in 0..ncols_q {
        let mut err2 = 0.0;
        let mut norm2 = 0.0;
        for i in 0..n {
            let want = (0..nrows_q).fold(T::ZERO, |acc, jr| acc + orig.col(1 + jr)[i] * q[c * ldq + jr]);
            err2 += (basis.col(1 + c)[i] - want).modulus_sq();
            norm2 += want.modulus_sq();
        }
        let rel = err2.sqrt() / norm2.sqrt().max(1e-300);
        assert!(rel < 1e-13, "Basis::mult_in_place col {c} rel err {rel:e}");
    }
}

#[test]
fn tile_mult_and_krylov_schur_agree_between_default_and_portable_kernels() {
    check_tile_mult_and_basis::<f64>();
    check_tile_mult_and_basis::<C64>();

    // Force Krylov-Schur restarts (`ncv = 12`, `nev = 4`, `block = 2`) so `k_tile_mult` executes
    // during thick-restart basis compression in both default (`faer`) and `portable-kernels` builds.
    // Shift `laplacian_1d(n) + 2*I` so eigenvalues `4 - 2 cos(k*pi/(n+1))` lie in `(2, 6)` and
    // relative residuals `||A x - lam x|| / (|lam| ||x||)` reach `< 1e-13`.
    let n = 60;
    let shift_diag = 2.0;
    let a = Csr::axpby(1.0, &laplacian_1d(n), shift_diag, &Csr::<f64>::identity(n));
    let b = Csr::<f64>::identity(n);
    let opts = EigenOptions { tol: 1e-15, block: 2, ncv: 12, par_solve_min_n: 0 };

    let real_pairs = eigs_near(
        &a,
        &b,
        1.9,
        4,
        Symmetry::Hermitian,
        Definiteness::PositiveDefinite,
        &opts,
        None,
    )
    .unwrap();
    assert!(real_pairs.converged);
    for (k, lam) in real_pairs.values.iter().enumerate() {
        let exact = shift_diag + 2.0 - 2.0 * ((k + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
        assert!((lam.re - exact).abs() / exact.max(1.0) < 1e-13, "mode {k}: {} vs {exact}", lam.re);
        assert!(real_pairs.residuals[k] < 1e-13, "eigs_near real residual {}: {:e}", k, real_pairs.residuals[k]);
    }

    let ac = a.to_complex();
    let bc = b.to_complex();
    let cplx_pairs = eigs_near(
        &ac,
        &bc,
        C64::new(1.9, 0.0),
        4,
        Symmetry::General,
        Definiteness::General,
        &opts,
        None,
    )
    .unwrap();
    assert!(cplx_pairs.converged);
    for (k, lam) in cplx_pairs.values.iter().enumerate() {
        let exact = shift_diag + 2.0 - 2.0 * ((k + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
        assert!((lam.re - exact).abs() / exact.max(1.0) < 1e-13, "cplx mode {k}: {lam} vs {exact}");
        assert!(cplx_pairs.residuals[k] < 1e-13, "eigs_near complex residual {}: {:e}", k, cplx_pairs.residuals[k]);
    }

    let rg_pairs = eigs_real_general_near(&a, &b, 1.9, 3, &opts).unwrap();
    assert!(rg_pairs.converged);
    for (k, &r) in rg_pairs.residuals.iter().enumerate() {
        assert!(r < 1e-13, "eigs_real_general_near residual {k}: {r:e}");
    }
}
