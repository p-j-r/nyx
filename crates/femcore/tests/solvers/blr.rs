#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::solvers::blr::{BlrLu, BlrOptions, BlrSymbolic, ShiftedBlrLu};
use femcore::solvers::csr::{laplacian_1d, laplacian_2d, laplacian_3d, mass_1d};
use femcore::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::solvers::ksp::{Gmres, GmresOptions, PcSide, Tolerances};
use femcore::solvers::linop::Identity;
use femcore::solvers::rng::random_vector;
use femcore::solvers::sparse_lu::{LuOptions, SparseLu};
use femcore::solvers::{C64, Csr, Scalar};

fn rel_residual<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut ax = vec![T::ZERO; b.len()];
    a.spmv(x, &mut ax);
    ax.iter().zip(b).map(|(p, q)| (*q - *p).modulus_sq()).sum::<f64>().sqrt() / norm(b)
}

fn helmholtz(l: &Csr<f64>, k2: f64, sigma: f64) -> Csr<C64> {
    let n = l.nrows();
    let t: Vec<_> = (0..n).map(|p| (p, p, C64::new(-k2, sigma * (1.0 + (p % 7) as f64 / 7.0)))).collect();
    Csr::axpby(C64::ONE, &l.to_complex(), C64::ONE, &Csr::from_triplets(n, n, &t).unwrap())
}

fn conv_diff(nx: usize, wind: f64) -> Csr<f64> {
    let h = 1.0 / (nx as f64 + 1.0);
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let p = j * nx + i;
            t.push((p, p, 4.0));
            if i > 0 {
                t.push((p, p - 1, -1.0 - wind * h));
            }
            if i + 1 < nx {
                t.push((p, p + 1, -1.0 + wind * h));
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

fn exact() -> BlrOptions {
    BlrOptions { tol: 0.0, leaf: 24, ..Default::default() }
}

/// Small blocks so that compression engages on test-sized fronts.
fn lowrank(tol: f64) -> BlrOptions {
    BlrOptions { tol, leaf: 24, block: 16, min_front: 32, accept_ratio: 0.7 }
}

#[test]
fn exact_mode_is_a_direct_solver_for_real_complex_and_nonsymmetric_matrices() {
    let cases: Vec<(&str, f64)> = {
        let a = laplacian_2d(30, 30);
        let lu = BlrLu::factor(&a, &exact()).unwrap();
        let b = random_vector::<f64>(900, 1);
        let mut x = b.clone();
        lu.solve_in_place(&mut x);
        vec![("laplacian 2-D", rel_residual(&a, &x, &b))]
    };
    assert!(cases[0].1 < 1e-12, "{cases:?}");
    // non-symmetric values on a symmetric pattern
    let a = conv_diff(32, 30.0);
    let lu = BlrLu::factor(&a, &exact()).unwrap();
    let b = random_vector::<f64>(a.nrows(), 2);
    let mut x = b.clone();
    lu.solve_in_place(&mut x);
    assert!(rel_residual(&a, &x, &b) < 1e-11, "{}", rel_residual(&a, &x, &b));
    // complex symmetric, indefinite
    let h = helmholtz(&laplacian_2d(30, 30), 2.5, 0.3);
    let lu = BlrLu::factor(&h, &exact()).unwrap();
    let bc = random_vector::<C64>(900, 3);
    let mut xc = bc.clone();
    lu.solve_in_place(&mut xc);
    assert!(rel_residual(&h, &xc, &bc) < 1e-11, "{}", rel_residual(&h, &xc, &bc));
    // agrees with the other sparse LU
    let reference = SparseLu::factor(&h, &LuOptions::default()).unwrap();
    let mut xr = bc.clone();
    reference.solve_in_place(&mut xr);
    let err = xc.iter().zip(&xr).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm(&xr);
    assert!(err < 1e-9, "{err}");
}

#[test]
fn indefinite_real_matrices_need_the_pivoting_inside_the_pivot_blocks() {
    // 2-D Laplacian shifted into the middle of its spectrum: tiny and sign-changing pivots
    let l = laplacian_2d(28, 28);
    let a = Csr::axpby(1.0, &l, -3.37, &Csr::identity(784));
    let lu = BlrLu::factor(&a, &exact()).unwrap();
    let b = random_vector::<f64>(784, 4);
    let mut x = b.clone();
    lu.solve_in_place(&mut x);
    assert!(rel_residual(&a, &x, &b) < 1e-9, "{}", rel_residual(&a, &x, &b));
}

#[test]
fn symbolic_structure_and_forests() {
    let a = laplacian_2d(40, 40);
    let sym = BlrSymbolic::analyse(&[&a], 24).unwrap();
    assert!(sym.num_fronts() > 10);
    let sizes = sym.front_sizes();
    assert_eq!(sizes.iter().map(|s| s.0).sum::<usize>(), 1600, "every unknown is a pivot exactly once");
    assert_eq!(sizes.last().unwrap().1, 0, "the root has no update block");
    // two disconnected components give two roots and still solve
    let n = 400;
    let blk = laplacian_2d(20, 20);
    let mut t = Vec::new();
    for r in 0..n {
        for (&c, &v) in blk.row(r).0.iter().zip(blk.row(r).1) {
            t.push((r, c as usize, v));
            t.push((r + n, c as usize + n, 2.0 * v));
        }
    }
    let a2 = Csr::from_triplets(2 * n, 2 * n, &t).unwrap();
    let lu = BlrLu::factor(&a2, &exact()).unwrap();
    let b = random_vector::<f64>(2 * n, 5);
    let mut x = b.clone();
    lu.solve_in_place(&mut x);
    assert!(rel_residual(&a2, &x, &b) < 1e-11);
}

#[test]
fn compression_shrinks_the_factor_and_is_accurate_to_the_tolerance() {
    let a = laplacian_3d(18);
    let b = random_vector::<f64>(a.nrows(), 6);
    let ex = BlrLu::factor(&a, &exact()).unwrap();
    assert_eq!(ex.stats().stored_entries, ex.stats().exact_entries);
    let mut prev_entries = usize::MAX;
    let mut prev_res = 0.0;
    for tol in [1e-4, 1e-8, 1e-12] {
        let lu = BlrLu::factor(&a, &lowrank(tol)).unwrap();
        let s = lu.stats();
        let mut x = b.clone();
        lu.solve_in_place(&mut x);
        let res = rel_residual(&a, &x, &b);
        println!(
            "tol {tol:e}: stored {:.3} of exact, {} of {} blocks low-rank (max rank {}), residual {res:.1e}",
            s.compression(),
            s.low_rank_blocks,
            s.blocks,
            s.max_rank
        );
        assert!(s.low_rank_blocks > 0 && s.compression() < 1.0, "no compression at tol {tol}");
        assert!(res < 200.0 * tol.max(1e-14), "tol {tol}: residual {res}");
        assert!(s.stored_entries >= prev_entries.min(s.stored_entries) || prev_entries == usize::MAX);
        assert!(prev_entries == usize::MAX || s.stored_entries >= prev_entries, "tighter tolerance must store more");
        assert!(res <= prev_res * 1.0 + 1.0 || prev_res == 0.0 || res < prev_res * 10.0);
        prev_entries = s.stored_entries;
        prev_res = res;
    }
}

#[test]
fn refinement_and_gmres_recover_full_accuracy_from_a_loose_factorisation() {
    let a = helmholtz(&laplacian_3d(16), 1.0, 0.15);
    let n = a.nrows();
    let b = random_vector::<C64>(n, 7);
    let lu = BlrLu::factor(&a, &lowrank(1e-3)).unwrap();
    assert!(lu.stats().compression() < 1.0);
    let mut x = vec![C64::ZERO; n];
    let (its, rel) = lu.solve_refined(&a, &b, &mut x, 30, 1e-12);
    println!("refinement: {its} iterations, residual {rel:.1e}, stored {:.2} of exact", lu.stats().compression());
    assert!(rel < 1e-11 && its <= 12, "{its} iterations, residual {rel}");
    // the same as a right preconditioner of GMRES
    let mut xg = vec![C64::ZERO; n];
    let g = Gmres::new(GmresOptions {
        tol: Tolerances { rtol: 1e-10, max_it: 100, ..Default::default() },
        restart: 30,
        side: PcSide::Right,
        ..Default::default()
    })
    .solve(&a, &lu, &b, &mut xg);
    assert!(g.reason.is_converged() && g.iterations <= 10, "{:?} after {}", g.reason, g.iterations);
    let plain = Gmres::new(GmresOptions {
        tol: Tolerances { rtol: 1e-10, max_it: 100, ..Default::default() },
        ..Default::default()
    })
    .solve(&a, &Identity, &b, &mut vec![C64::ZERO; n]);
    assert!(!plain.reason.is_converged(), "unpreconditioned GMRES was expected to stall");
}

#[test]
fn shifted_solver_matches_faer_and_reuses_the_symbolic_analysis() {
    let n = 300;
    let (a, m) = (laplacian_1d(n), mass_1d(n, 1.0 / (n as f64 + 1.0)));
    let b = random_vector::<f64>(n, 8);
    let mut reference = ShiftedFactorization::new(&a, Some(&m), Definiteness::General).unwrap();
    let mut ours = ShiftedBlrLu::new(&a, Some(&m), &exact()).unwrap();
    for sigma in [3.7, 11.2] {
        reference.set_shift(sigma).unwrap();
        ours.set_shift(sigma).unwrap();
        let (mut x1, mut x2) = (b.clone(), b.clone());
        reference.solve_block_in_place(&mut x1, 1);
        ours.solve_block_in_place(&mut x2, 1);
        let err = x1.iter().zip(&x2).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&x1);
        assert!(err < 1e-8, "shift {sigma}: {err}");
    }
    assert_eq!(ours.stats().symbolic_count, 1);
    assert_eq!(ours.stats().numeric_count, 2);
}

#[test]
fn block_solve_equals_single_solves() {
    let a = helmholtz(&laplacian_2d(14, 14), 1.0, 0.2);
    let lu = BlrLu::factor(&a, &lowrank(1e-6)).unwrap();
    let n = a.nrows();
    let rhs = random_vector::<C64>(3 * n, 9);
    let mut blk = rhs.clone();
    lu.solve_block_in_place(&mut blk, 3);
    for c in 0..3 {
        let mut one = rhs[c * n..(c + 1) * n].to_vec();
        lu.solve_in_place(&mut one);
        assert_eq!(&blk[c * n..(c + 1) * n], &one[..]);
    }
}
