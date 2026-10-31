//! S4 study: the block low-rank multifrontal LU (`blr.rs`) on 2-D and 3-D grid operators: how much of the exact
//! multifrontal factor is stored, how fast factorisation and solve are, and how accurate the compressed solver is
//! (direct, after iterative refinement, and as a GMRES preconditioner). `cargo run --release --example blr_study`
//! Optional arguments: `n2d n3d block min_front`.
use femcore::solvers::blr::{BlrLu, BlrOptions};
use femcore::solvers::csr::{laplacian_2d, laplacian_3d};
use femcore::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::solvers::ksp::{Gmres, GmresOptions, PcSide, Tolerances};
use femcore::solvers::rng::random_vector;
use femcore::solvers::{C64, Csr};
use std::time::Instant;

fn helmholtz(l: &Csr<f64>, k2: f64, sigma: f64) -> Csr<C64> {
    let n = l.nrows();
    let t: Vec<_> = (0..n).map(|p| (p, p, C64::new(-k2, sigma * (1.0 + (p % 7) as f64 / 7.0)))).collect();
    Csr::axpby(C64::ONE, &l.to_complex(), C64::ONE, &Csr::from_triplets(n, n, &t).unwrap())
}

fn residual(a: &Csr<C64>, x: &[C64], b: &[C64]) -> f64 {
    let mut ax = vec![C64::ZERO; b.len()];
    a.spmv(x, &mut ax);
    ax.iter().zip(b).map(|(p, q)| (*q - *p).norm_sqr()).sum::<f64>().sqrt()
        / b.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt()
}

fn study(name: &str, a: &Csr<C64>, block: usize, min_front: usize) {
    let n = a.nrows();
    let b = random_vector::<C64>(n, 1);
    println!("\n{name}: n = {n}");
    println!(
        "{:>8} {:>8} {:>10} {:>10} {:>9} {:>9} {:>10} {:>8} {:>10} {:>9}",
        "tol", "stored", "factor[s]", "solve[ms]", "low/all", "maxrank", "residual", "refine", "gmres its", "largest"
    );
    for tol in [0.0, 1e-4, 1e-8] {
        let o = BlrOptions { tol, leaf: 64, block, min_front, accept_ratio: 0.7 };
        let t = Instant::now();
        let lu = BlrLu::factor(a, &o).unwrap();
        let tf = t.elapsed().as_secs_f64();
        let mut x = b.clone();
        let t = Instant::now();
        lu.solve_in_place(&mut x);
        let ts = t.elapsed().as_secs_f64() * 1e3;
        let res = residual(a, &x, &b);
        let mut xr = vec![C64::ZERO; n];
        let (its, _) = lu.solve_refined(a, &b, &mut xr, 30, 1e-12);
        let mut xg = vec![C64::ZERO; n];
        let g = Gmres::new(GmresOptions {
            tol: Tolerances { rtol: 1e-10, max_it: 200, ..Default::default() },
            side: PcSide::Right,
            restart: 40,
            ..Default::default()
        })
        .solve(a, &lu, &b, &mut xg);
        let s = lu.stats();
        println!(
            "{tol:>8.0e} {:>8.3} {tf:>10.2} {ts:>10.2} {:>4}/{:<4} {:>9} {res:>10.1e} {its:>8} {:>10} {:>9}",
            s.compression(),
            s.low_rank_blocks,
            s.blocks,
            s.max_rank,
            g.iterations,
            format!("{}x{}", s.largest_front.0, s.largest_front.1)
        );
    }
    let mut f = ShiftedFactorization::new(a, None, Definiteness::General).unwrap();
    f.set_shift(C64::ZERO).unwrap();
    let st = f.stats();
    let mut x = b.clone();
    let t = Instant::now();
    f.solve_block_in_place(&mut x, 1);
    println!(
        "faer LU (reference): factor {:.2} s, solve {:.2} ms, residual {:.1e}",
        st.numeric_secs_last,
        t.elapsed().as_secs_f64() * 1e3,
        residual(a, &x, &b)
    );
}

fn main() {
    let args: Vec<usize> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let (n2, n3) = (args.first().copied().unwrap_or(300), args.get(1).copied().unwrap_or(30));
    let (block, min_front) = (args.get(2).copied().unwrap_or(64), args.get(3).copied().unwrap_or(128));
    study("2-D 5-point Helmholtz + absorption", &helmholtz(&laplacian_2d(n2, n2), 1.0, 0.1), block, min_front);
    study("3-D 7-point Helmholtz + absorption", &helmholtz(&laplacian_3d(n3), 1.0, 0.1), block, min_front);
}
