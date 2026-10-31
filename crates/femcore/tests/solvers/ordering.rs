//! Fill-reducing orderings for the sparse Cholesky: all choices give the same solution; nested dissection
//! reduces fill on 3-D graphs; `Auto` keeps the better of the two.
use crate::common;

use common::*;
use femcore::solvers::Csr;
use femcore::solvers::csr::{laplacian_2d, laplacian_3d};
use femcore::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::solvers::ordering::FillOrdering;
use femcore::solvers::rng::random_vector;
use std::time::Instant;

fn factor(a: &Csr<f64>, ord: FillOrdering) -> (ShiftedFactorization<f64>, f64) {
    let mut f = ShiftedFactorization::new(a, None, Definiteness::PositiveDefinite).unwrap();
    f.ordering = ord;
    let t = Instant::now();
    f.set_shift(0.0).unwrap();
    (f, t.elapsed().as_secs_f64())
}

fn solve(f: &ShiftedFactorization<f64>, b: &[f64]) -> Vec<f64> {
    let mut x = b.to_vec();
    f.solve_block_in_place(&mut x, 1);
    x
}

#[test]
fn every_ordering_gives_the_same_solution() {
    for a in [laplacian_3d(16), laplacian_2d(70, 70)] {
        let b = random_vector::<f64>(a.nrows(), 1);
        let reference = solve(&factor(&a, FillOrdering::Amd).0, &b);
        for ord in [FillOrdering::NestedDissection, FillOrdering::Auto] {
            let x = solve(&factor(&a, ord).0, &b);
            let err = x.iter().zip(&reference).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&reference);
            assert!(err < 1e-11, "{ord:?}: {err}");
        }
        // residual
        let mut r = vec![0.0; a.nrows()];
        a.spmv(&reference, &mut r);
        assert!(r.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-10));
    }
}

#[test]
fn nested_dissection_reduces_fill_on_3d_graphs_and_auto_picks_the_smaller() {
    let a3 = laplacian_3d(32);
    let (amd, nd, auto) = (
        factor(&a3, FillOrdering::Amd).0,
        factor(&a3, FillOrdering::NestedDissection).0,
        factor(&a3, FillOrdering::Auto).0,
    );
    let (na, nn, nu) = (amd.stats().factor_nnz, nd.stats().factor_nnz, auto.stats().factor_nnz);
    eprintln!("3-D 32³: factor nnz AMD {na}, ND {nn}, Auto {nu}");
    assert!(nn < na, "ND should beat AMD on a 3-D grid: {nn} vs {na}");
    assert_eq!(nu, nn.min(na));
    // 2-D: AMD has less fill, Auto must not choose worse than AMD
    let a2 = laplacian_2d(120, 120);
    let (amd2, nd2, auto2) = (
        factor(&a2, FillOrdering::Amd).0,
        factor(&a2, FillOrdering::NestedDissection).0,
        factor(&a2, FillOrdering::Auto).0,
    );
    eprintln!(
        "2-D 120²: factor nnz AMD {}, ND {}, Auto {}",
        amd2.stats().factor_nnz,
        nd2.stats().factor_nnz,
        auto2.stats().factor_nnz
    );
    assert!(auto2.stats().factor_nnz <= amd2.stats().factor_nnz.min(nd2.stats().factor_nnz));
}

#[test]
#[ignore = "timing study; run with --ignored --nocapture"]
fn timing_3d() {
    let a = laplacian_3d(40);
    let b = random_vector::<f64>(a.nrows(), 1);
    use femcore::solvers::direct::Threads;
    for (ord, th) in [FillOrdering::Amd, FillOrdering::NestedDissection, FillOrdering::Auto]
        .into_iter()
        .flat_map(|o| [(o, Threads::Global), (o, Threads::Sequential)])
    {
        let mut f = ShiftedFactorization::new(&a, None, Definiteness::PositiveDefinite).unwrap();
        f.ordering = ord;
        f.threads.factor = th;
        let t = Instant::now();
        f.set_shift(0.0).unwrap();
        let t = t.elapsed().as_secs_f64();
        eprintln!("{th:?} factor threads:");
        let t0 = Instant::now();
        let _ = solve(&f, &b);
        eprintln!(
            "{ord:?}: nnz(L) {}, set_shift {t:.2} s (analysis {:.2} s, numeric {:.2} s), one solve {:.3} s",
            f.stats().factor_nnz,
            f.stats().analyse_secs,
            f.stats().numeric_secs_last,
            t0.elapsed().as_secs_f64()
        );
    }
}
