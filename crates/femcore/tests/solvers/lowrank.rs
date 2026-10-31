#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::solvers::csr::{laplacian_1d, laplacian_2d};
use femcore::solvers::ksp::{Gmres, GmresOptions, Tolerances};
use femcore::solvers::linop::Identity;
use femcore::solvers::lowrank::{Hodlr, HodlrOptions, aca_full, aca_partial};
use femcore::solvers::rng::random_vector;
use femcore::solvers::sparse_lu::{LuOptions, SparseLu};
use femcore::solvers::{C64, Csr, Scalar};

fn rel_err<T: Scalar>(a: &[T], b: &[T]) -> f64 {
    a.iter().zip(b).map(|(p, q)| (*p - *q).modulus_sq()).sum::<f64>().sqrt() / norm(b)
}

/// 1-D log kernel with a dominant diagonal: SPD-ish, well conditioned, off-diagonal blocks of low rank.
fn log_kernel(n: usize) -> impl Fn(usize, usize) -> f64 {
    move |i, j| {
        if i == j {
            2.0 + n as f64 * 0.01
        } else {
            let d = (i as f64 - j as f64).abs() / n as f64;
            d.ln() * 0.05
        }
    }
}

fn dense_apply<T: Scalar>(entry: &dyn Fn(usize, usize) -> T, n: usize, x: &[T]) -> Vec<T> {
    (0..n).map(|i| (0..n).fold(T::ZERO, |s, j| s + entry(i, j) * x[j])).collect()
}

#[test]
fn aca_recovers_exact_low_rank_blocks_real_and_complex() {
    let (m, n, r) = (60, 45, 4);
    let uf = random_vector::<f64>(m * r, 1);
    let vf = random_vector::<f64>(n * r, 2);
    let entry = |i: usize, j: usize| (0..r).map(|k| uf[k * m + i] * vf[k * n + j]).sum::<f64>();
    let lr = aca_partial(m, n, &entry, 1e-12, 30);
    assert!(lr.r <= r + 2, "rank {}", lr.r);
    let d = lr.to_dense();
    let mut err = 0.0f64;
    for j in 0..n {
        for i in 0..m {
            err = err.max((d[j * m + i] - entry(i, j)).abs());
        }
    }
    assert!(err < 1e-10, "{err}");
    // complex
    let uc = random_vector::<C64>(m * r, 3);
    let vc = random_vector::<C64>(n * r, 4);
    let entry_c = |i: usize, j: usize| (0..r).fold(C64::ZERO, |s, k| s + uc[k * m + i] * vc[k * n + j]);
    let lr = aca_partial(m, n, &entry_c, 1e-12, 30);
    assert!(lr.r <= r + 2, "complex rank {}", lr.r);
    let d = lr.to_dense();
    let mut err = 0.0f64;
    for j in 0..n {
        for i in 0..m {
            err = err.max((d[j * m + i] - entry_c(i, j)).norm());
        }
    }
    assert!(err < 1e-10, "{err}");
    // full pivoting on the same data
    let mut dense: Vec<C64> = (0..m * n).map(|k| entry_c(k % m, k / m)).collect();
    let full = aca_full(&mut dense, m, n, 1e-13, 30);
    assert_eq!(full.r, r);
    // a zero block has rank 0
    assert_eq!(aca_partial(10, 10, &|_, _| 0.0f64, 1e-8, 5).r, 0);
}

#[test]
fn rank_grows_slowly_as_the_tolerance_shrinks() {
    let n = 400;
    let k = log_kernel(n);
    let block = |i: usize, j: usize| k(i, 200 + j); // far off-diagonal block 200×200
    let ranks: Vec<usize> =
        [1e-3, 1e-6, 1e-9, 1e-12].iter().map(|&t| aca_partial(200, 200, &block, t, 100).r).collect();
    println!("ranks for tol 1e-3/1e-6/1e-9/1e-12: {ranks:?}");
    assert!(ranks.windows(2).all(|w| w[0] <= w[1]));
    assert!(ranks[3] < 40, "{ranks:?}");
}

#[test]
fn hodlr_matvec_and_solve_on_a_kernel_matrix() {
    let n = 1024;
    let k = log_kernel(n);
    let tol = 1e-8;
    let h = Hodlr::from_entries(n, &k, HodlrOptions { leaf: 32, tol, max_rank: 100 });
    println!("rank profile {:?}, storage {:.3} of dense", h.rank_profile(), h.storage() as f64 / (n * n) as f64);
    assert!(h.storage() * 3 < n * n);
    let x = random_vector::<f64>(n, 5);
    let exact = dense_apply(&k, n, &x);
    let mut y = vec![0.0; n];
    h.matvec(&x, &mut y);
    assert!(rel_err(&y, &exact) < 100.0 * tol, "matvec error {}", rel_err(&y, &exact));
    let lu = h.factor().expect("factorisation");
    let b = random_vector::<f64>(n, 6);
    let mut sol = b.clone();
    lu.solve_in_place(&mut sol);
    let back = dense_apply(&k, n, &sol);
    assert!(rel_err(&back, &b) < 1000.0 * tol, "solve residual {}", rel_err(&back, &b));
}

#[test]
fn hodlr_complex_kernel() {
    let n = 600;
    let entry = move |i: usize, j: usize| -> C64 {
        if i == j {
            C64::new(3.0, 0.5)
        } else {
            let d = (i as f64 - j as f64).abs() / n as f64;
            C64::from_polar(0.05 / (0.1 + d), 6.0 * d)
        }
    };
    let h = Hodlr::from_entries(n, &entry, HodlrOptions { leaf: 24, tol: 1e-9, max_rank: 120 });
    let lu = h.factor().unwrap();
    let b = random_vector::<C64>(n, 7);
    let mut x = b.clone();
    lu.solve_in_place(&mut x);
    let back = dense_apply(&entry, n, &x);
    assert!(rel_err(&back, &b) < 1e-6, "{}", rel_err(&back, &b));
}

#[test]
fn sparse_banded_matrices_are_exactly_hodlr_with_rank_one() {
    let a = laplacian_1d(500);
    let h = Hodlr::from_csr(&a, HodlrOptions { leaf: 20, tol: 1e-13, max_rank: 10 });
    assert!(h.rank_profile().iter().all(|&r| r == 1), "{:?}", h.rank_profile());
    let b = random_vector::<f64>(500, 8);
    let lu = h.factor().unwrap();
    let mut x = b.clone();
    lu.solve_in_place(&mut x);
    let mut ax = vec![0.0; 500];
    a.spmv(&x, &mut ax);
    assert!(rel_err(&ax, &b) < 1e-9, "{}", rel_err(&ax, &b));
}

/// The motivating case: the Schur complement of a 2-D grid onto a separator line is dense, and its
/// off-diagonal blocks have low numerical rank.
#[test]
fn schur_complement_on_a_separator_is_hodlr_compressible() {
    let nx = 90;
    let a = laplacian_2d(nx, nx);
    let sep_col = nx / 2;
    let is_sep = |p: usize| p % nx == sep_col;
    let interior: Vec<usize> = (0..nx * nx).filter(|&p| !is_sep(p)).collect();
    let sep: Vec<usize> = (0..nx * nx).filter(|&p| is_sep(p)).collect();
    let mut pos = vec![usize::MAX; nx * nx];
    for (k, &p) in interior.iter().enumerate() {
        pos[p] = k;
    }
    let mut t_ii = Vec::new();
    let mut a_is: Vec<Vec<(usize, f64)>> = vec![Vec::new(); sep.len()];
    let mut sep_pos = vec![usize::MAX; nx * nx];
    for (k, &p) in sep.iter().enumerate() {
        sep_pos[p] = k;
    }
    for &p in &interior {
        let (cols, vals) = a.row(p);
        for (&c, &v) in cols.iter().zip(vals) {
            let c = c as usize;
            if is_sep(c) {
                a_is[sep_pos[c]].push((pos[p], v));
            } else {
                t_ii.push((pos[p], pos[c], v));
            }
        }
    }
    let nint = interior.len();
    let a_ii = Csr::from_triplets(nint, nint, &t_ii).unwrap();
    let lu = SparseLu::factor(&a_ii, &LuOptions::default()).unwrap();
    let ns = sep.len();
    let mut s = vec![0.0f64; ns * ns]; // column-major Schur complement
    for (jj, &pj) in sep.iter().enumerate() {
        let mut col = vec![0.0; nint];
        for &(i, v) in &a_is[jj] {
            col[i] = v;
        }
        lu.solve_in_place(&mut col); // A_ii⁻¹ A_is e_j
        // S[:, j] = A_ss[:, j] − A_si A_ii⁻¹ A_is e_j
        for (ii, &pi) in sep.iter().enumerate() {
            let mut v = a.row(pi).0.iter().zip(a.row(pi).1).find(|(c, _)| **c as usize == pj).map_or(0.0, |(_, v)| *v);
            for &(i, w) in &a_is[ii] {
                v -= w * col[i];
            }
            s[jj * ns + ii] = v;
        }
    }
    let h = Hodlr::from_dense(&s, ns, HodlrOptions { leaf: 8, tol: 1e-8, max_rank: ns });
    let profile = h.rank_profile();
    println!(
        "Schur complement {ns}×{ns}: ranks per level {profile:?}, storage {:.2} of dense",
        h.storage() as f64 / (ns * ns) as f64
    );
    assert!(profile[0] <= ns / 4, "top-level rank {} of block size {}", profile[0], ns / 2);
    // and it is a faithful representation
    let x = random_vector::<f64>(ns, 9);
    let mut y = vec![0.0; ns];
    h.matvec(&x, &mut y);
    let mut exact = vec![0.0; ns];
    for j in 0..ns {
        for i in 0..ns {
            exact[i] += s[j * ns + i] * x[j];
        }
    }
    assert!(rel_err(&y, &exact) < 1e-6);
}

#[test]
fn a_loose_hodlr_factorisation_is_an_effective_preconditioner() {
    let n = 800;
    let entry = move |i: usize, j: usize| -> f64 {
        if i == j {
            1.0
        } else {
            let d = (i as f64 - j as f64).abs() / n as f64;
            0.4 / (0.02 + d) / n as f64 * 5.0
        }
    };
    let mut dense = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            dense[j * n + i] = entry(i, j);
        }
    }
    let mut trip = Vec::new();
    for i in 0..n {
        for j in 0..n {
            trip.push((i, j, dense[j * n + i]));
        }
    }
    let a = Csr::from_triplets(n, n, &trip).unwrap();
    let b = random_vector::<f64>(n, 10);
    let tol = Tolerances { rtol: 1e-10, max_it: 2000, ..Default::default() };
    let mut x0 = vec![0.0; n];
    let plain = Gmres::new(GmresOptions { tol, restart: 60, ..Default::default() }).solve(&a, &Identity, &b, &mut x0);
    let h = Hodlr::from_dense(&dense, n, HodlrOptions { leaf: 32, tol: 1e-2, max_rank: 20 });
    let pre = h.factor().unwrap();
    let mut x1 = vec![0.0; n];
    let prec = Gmres::new(GmresOptions { tol, restart: 60, ..Default::default() }).solve(&a, &pre, &b, &mut x1);
    println!(
        "GMRES(60): {} its plain, {} its with a tol 1e-2 HODLR preconditioner (storage {:.2})",
        plain.iterations,
        prec.iterations,
        h.storage() as f64 / (n * n) as f64
    );
    assert!(prec.reason.is_converged());
    assert!(prec.iterations <= 12, "{}", prec.iterations);
    assert!(prec.iterations * 3 <= plain.iterations.max(36), "plain {}", plain.iterations);
}
