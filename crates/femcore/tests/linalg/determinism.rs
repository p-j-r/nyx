#![allow(clippy::needless_range_loop)] // index loops mirror the mathematics in these tests
//! Results must not depend on the number of threads (fixed chunking, ordered folds).


use femcore::linalg::csr::{laplacian_2d, mass_1d};
use femcore::linalg::eps::{EpsOptions, krylov_schur};
use femcore::linalg::{C64, Csr, LinearOperator, Scalar, vecops};

fn in_pool<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(f)
}

fn run<T: Scalar + femcore::linalg::faer_glue::FaerScalar>(
    a: &Csr<T>,
    weight: Option<&dyn LinearOperator<T>>,
    opts: &EpsOptions,
) -> (Vec<femcore::linalg::C64>, Vec<T>, usize) {
    let r = krylov_schur(a, weight, None, opts).unwrap();
    let vecs: Vec<T> = (0..r.values.len()).flat_map(|i| r.vector(i).to_vec()).collect();
    (r.values, vecs, r.matvecs)
}

#[test]
fn krylov_schur_is_bitwise_independent_of_thread_count() {
    // n is large enough that the parallel thresholds engage.
    let a = laplacian_2d(160, 150);
    let opts = EpsOptions { tol: 1e-4, ncv: Some(24), max_it: Some(60), block_size: 2, ..EpsOptions::hermitian(4) };
    let one = in_pool(1, || run(&a, None, &opts));
    let four = in_pool(4, || run(&a, None, &opts));
    assert_eq!(one.2, four.2);
    assert_eq!(one.0, four.0, "eigenvalues differ between thread counts");
    assert!(one.1 == four.1, "eigenvectors differ between thread counts");
}

#[test]
fn weighted_complex_general_is_bitwise_independent_of_thread_count() {
    let n = 30_000;
    let a = Csr::axpby(C64::ONE, &common_lap(n), C64::new(0.0, 0.1), &Csr::identity(n));
    let opts = EpsOptions { tol: 1e-6, ncv: Some(20), max_it: Some(60), ..EpsOptions::general(3) };
    let one = in_pool(1, || run(&a, None, &opts));
    let four = in_pool(4, || run(&a, None, &opts));
    assert_eq!(one.0, four.0);
    assert!(one.1 == four.1);
    // Weighted Hermitian path: Op = M is trivially self-adjoint in the M inner product, which
    // keeps the test independent of a factorisation. (Trait objects are not `Send`, so each
    // run builds its own.)
    let opts = EpsOptions { tol: 1e-6, ncv: Some(20), max_it: Some(60), ..EpsOptions::hermitian(3) };
    let weighted = || {
        let m = mass_1d(n, 1.0).to_complex();
        let w: &dyn LinearOperator<C64> = &m;
        run(&m, Some(w), &opts)
    };
    let one = in_pool(1, weighted);
    let four = in_pool(4, weighted);
    assert_eq!(one.0, four.0);
    assert!(one.1 == four.1);
}

fn common_lap(n: usize) -> Csr<C64> {
    femcore::linalg::csr::laplacian_1d(n).to_complex()
}

#[test]
fn reductions_are_bitwise_independent_of_thread_count() {
    let n = 300_001;
    let x = femcore::linalg::rng::random_vector::<C64>(n, 1);
    let y = femcore::linalg::rng::random_vector::<C64>(n, 2);
    let d1 = in_pool(1, || vecops::dotc(&x, &y));
    let d4 = in_pool(4, || vecops::dotc(&x, &y));
    assert_eq!(d1, d4);
    let a = laplacian_2d(400, 300).to_complex();
    let x = femcore::linalg::rng::random_vector::<C64>(a.nrows(), 3);
    let (mut y1, mut y4) = (vec![C64::ZERO; a.nrows()], vec![C64::ZERO; a.nrows()]);
    let s1 = in_pool(1, || a.spmv_dotc(&x, &mut y1));
    let s4 = in_pool(4, || a.spmv_dotc(&x, &mut y4));
    assert_eq!(s1, s4);
    assert_eq!(y1, y4);
}

#[test]
fn matrix_free_operator_is_bitwise_independent_of_thread_count() {
    use femcore::linalg::matfree::{Boundary, SpectralBox};
    // 14³ elements of degree 4: 57³ = 185 193 nodes, above the parallel threshold
    let op = SpectralBox::<f64>::uniform(4, &[14, 14, 14], &[0.1, 0.1, 0.1], 1.0, -3.0, Boundary::Dirichlet);
    let x = femcore::linalg::rng::random_vector::<f64>(op.len(), 3);
    let run = || {
        let mut y = vec![0.0; op.len()];
        op.apply(&x, &mut y);
        y
    };
    let one = in_pool(1, run);
    let four = in_pool(4, run);
    assert!(one == four, "matrix-free product differs between thread counts");
}
