//! S2/LU experiment: what does a user-chosen fill-reducing ordering buy for *LU* (the matrices faer's LU
//! orders with its own COLAMD)? Compares `femcore::linalg::sparse_lu::SparseLu` with natural / AMD / nested
//! dissection against faer's supernodal sparse LU (`ShiftedFactorization`, `Definiteness::General`) on a 2-D
//! and a 3-D complex-symmetric Helmholtz-type matrix: factor entries, factorisation and solve time, residual.
//! `cargo run --release --example lu_ordering_study [n2d n3d]`
use femcore::linalg::csr::{laplacian_2d, laplacian_3d};
use femcore::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::linalg::rng::random_vector;
use femcore::linalg::sparse_lu::{LuOptions, LuOrdering, SparseLu};
use femcore::linalg::{C64, Csr};
use std::time::Instant;

fn shifted(l: &Csr<f64>, k2: f64, sigma: f64) -> Csr<C64> {
    let n = l.nrows();
    let d: Vec<_> = (0..n).map(|p| (p, p, C64::new(-k2, sigma * (1.0 + (p % 7) as f64 / 7.0)))).collect();
    Csr::axpby(C64::ONE, &l.to_complex(), C64::ONE, &Csr::from_triplets(n, n, &d).unwrap())
}

fn residual(a: &Csr<C64>, x: &[C64], b: &[C64]) -> f64 {
    let mut ax = vec![C64::ZERO; b.len()];
    a.spmv(x, &mut ax);
    let r: f64 = ax.iter().zip(b).map(|(p, q)| (*q - *p).norm_sqr()).sum::<f64>().sqrt();
    r / b.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt()
}

fn study(name: &str, a: &Csr<C64>) {
    let n = a.nrows();
    let b = random_vector::<C64>(n, 1);
    println!("\n{name}: n = {n}, nnz(A) = {}", a.nnz());
    println!(
        "{:<22} {:>12} {:>10} {:>10} {:>10} {:>10}",
        "method", "nnz(L+U)", "order [s]", "factor [s]", "solve [s]", "residual"
    );
    for (label, ordering) in [
        ("sparse_lu natural", LuOrdering::Natural),
        ("sparse_lu amd", LuOrdering::Amd),
        ("sparse_lu nd", LuOrdering::NestedDissection { leaf: 64 }),
    ] {
        let t = Instant::now();
        let lu = SparseLu::factor(a, &LuOptions { ordering, pivot_tol: 0.1 }).unwrap();
        let total = t.elapsed().as_secs_f64();
        let mut x = b.clone();
        let ts = Instant::now();
        lu.solve_in_place(&mut x);
        println!(
            "{label:<22} {:>12} {:>10.3} {:>10.3} {:>10.4} {:>10.1e}",
            lu.factor_nnz(),
            lu.stats().order_secs,
            total - lu.stats().order_secs,
            ts.elapsed().as_secs_f64(),
            residual(a, &x, &b)
        );
    }
    let mut f = ShiftedFactorization::new(a, None, Definiteness::General).unwrap();
    f.set_shift(C64::ZERO).unwrap();
    let st = f.stats();
    let mut x = b.clone();
    let ts = Instant::now();
    f.solve_block_in_place(&mut x, 1);
    println!(
        "{:<22} {:>12} {:>10.3} {:>10.3} {:>10.4} {:>10.1e}",
        "faer LU (colamd)",
        "n/a",
        st.analyse_secs,
        st.numeric_secs_last,
        ts.elapsed().as_secs_f64(),
        residual(a, &x, &b)
    );
}

fn main() {
    let args: Vec<usize> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let (n2, n3) = (args.first().copied().unwrap_or(160), args.get(1).copied().unwrap_or(24));
    study("2-D 5-point Helmholtz + absorption", &shifted(&laplacian_2d(n2, n2), 1.5, 0.3));
    study("3-D 7-point Helmholtz + absorption", &shifted(&laplacian_3d(n3), 1.5, 0.3));
}
