// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Cost per right-hand side of a sparse triangular solve as a function of the block width.
//! The solves dominate shift-and-invert eigensolves, so this decides whether the block variant
//! of Krylov–Schur is also *faster* per vector, not just more robust for degenerate spectra.
//!
//!   cargo run --release --example solve_scaling -- [n2d] [n3d]

use femcore::linalg::csr::{Csr, laplacian_2d, laplacian_3d};
use femcore::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, Threads};
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Scalar};
use std::time::Instant;

fn run<T: femcore::linalg::faer_glue::FaerScalar>(label: &str, a: &Csr<T>, def: Definiteness, sigma: T) {
    let n = a.nrows();
    let mut f = ShiftedFactorization::new(a, None, def).unwrap();
    // FEMSOLVE_SOLVE_THREADS=global|seq|<n> overrides the (sequential) default solve policy.
    f.threads.solve = match std::env::var("FEMSOLVE_SOLVE_THREADS").as_deref() {
        Ok("global") => Threads::Global,
        Ok("seq") | Err(_) => Threads::Sequential,
        Ok(n) => Threads::Fixed(n.parse().expect("number of threads")),
    };
    let t0 = Instant::now();
    f.set_shift(sigma).unwrap();
    println!(
        "{label}: n = {n}, factorisation {:.3}s ({:?}), solve threads {:?}",
        t0.elapsed().as_secs_f64(),
        f.stats().used,
        f.threads.solve
    );
    let mut base = 0.0;
    for &nrhs in &[1usize, 2, 4, 8] {
        let mut rhs = random_vector::<T>(n * nrhs, 1);
        f.solve_block_in_place(&mut rhs, nrhs); // warm-up
        let mut best = f64::MAX;
        for _ in 0..5 {
            let mut r = random_vector::<T>(n * nrhs, 2);
            let t = Instant::now();
            f.solve_block_in_place(&mut r, nrhs);
            best = best.min(t.elapsed().as_secs_f64());
        }
        let per = best / nrhs as f64;
        if nrhs == 1 {
            base = per;
        }
        println!(
            "  nrhs {nrhs}: {:8.3} ms total, {:8.3} ms per vector ({:.2}x of single)",
            best * 1e3,
            per * 1e3,
            per / base
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n2d: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(200);
    let n3d: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);
    let a2 = Csr::axpby(C64::ONE, &laplacian_2d(n2d, n2d).to_complex(), C64::new(0.0, 0.3), &Csr::identity(n2d * n2d));
    run("2-D complex symmetric, LU", &a2, Definiteness::General, C64::new(0.5, 0.01));
    run("3-D real SPD, Cholesky", &laplacian_3d(n3d), Definiteness::PositiveDefinite, -0.1);
    run("3-D real SPD, LU", &laplacian_3d(n3d), Definiteness::General, -0.1);
    let _ = <f64 as Scalar>::ONE;
}
