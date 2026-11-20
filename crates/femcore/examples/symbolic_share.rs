// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! How much of one factorisation is *symbolic* (ordering and elimination structure), i.e. how much
//! symbolic reuse (port item S1) can save per additional shift.
//!
//!   cargo run --release --example symbolic_share -- [n2d] [n3d]
//!
//! Matrices follow `examples/lu_scaling.rs` of nyx as described in the port list: a complex
//! symmetric 2-D problem (LU) and a real SPD 3-D problem (Cholesky and, for comparison, LU).

use femcore::linalg::csr::{Csr, laplacian_2d, laplacian_3d};
use femcore::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::linalg::{C64, Scalar};
use std::time::Instant;

fn report<T: femcore::linalg::faer_glue::FaerScalar>(label: &str, a: &Csr<T>, def: Definiteness, shifts: &[T]) {
    let t0 = Instant::now();
    let mut f = ShiftedFactorization::new(a, None, def).expect("setup");
    let setup = t0.elapsed().as_secs_f64();
    let mut first = 0.0;
    let mut later = Vec::new();
    for (i, &s) in shifts.iter().enumerate() {
        let t = Instant::now();
        f.set_shift(s).expect("factorise");
        let dt = t.elapsed().as_secs_f64();
        if i == 0 {
            first = dt;
        } else {
            later.push(dt);
        }
    }
    let st = f.stats();
    let later_mean = later.iter().sum::<f64>() / later.len().max(1) as f64;
    let symbolic = st.analyse_secs - setup; // analyses done lazily inside the first set_shift
    println!(
        "{label:<34} n={:>7} nnz={:>8} | pattern {setup:>6.3}s symbolic {symbolic:>6.3}s | 1st shift {first:>7.3}s, later shifts {later_mean:>7.3}s (numeric {:>7.3}s) | symbolic share of 1st: {:>5.1}% | factor {:?}",
        a.nrows(),
        a.nnz(),
        st.numeric_secs_last,
        100.0 * symbolic / first.max(1e-12),
        st.used
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n2d: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(250);
    let n3d: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);

    let a2 = Csr::axpby(C64::ONE, &laplacian_2d(n2d, n2d).to_complex(), C64::new(0.0, 0.3), &Csr::identity(n2d * n2d));
    let shifts_c: Vec<C64> = (0..4).map(|i| C64::new(0.5 + 0.1 * i as f64, 0.01)).collect();
    report("2-D complex symmetric, LU", &a2, Definiteness::General, &shifts_c);

    let a3 = laplacian_3d(n3d);
    let shifts_r: Vec<f64> = (0..4).map(|i| -0.1 - 0.05 * i as f64).collect();
    report("3-D real SPD, Cholesky", &a3, Definiteness::PositiveDefinite, &shifts_r);
    report("3-D real SPD, LU (same matrix)", &a3, Definiteness::General, &shifts_r);
    let _ = <f64 as Scalar>::ONE;
}
