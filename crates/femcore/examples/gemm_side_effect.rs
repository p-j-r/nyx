// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Does one call into faer's GEMM slow down *later, unrelated* code on this CPU?
//! (Observed on the benchmark VM: after the first restart multiply every phase of an eigensolve
//! ran about 3x slower, including code that never calls faer. A dirty upper vector state after
//! AVX-512 kernels, or a clock/licence effect, would explain it.)
//!
//!   cargo run --release --example gemm_side_effect
//!   cargo run --release --example gemm_side_effect --features portable-kernels   # control

use femcore::linalg::C64;
use femcore::linalg::basis::Basis;
use femcore::linalg::rng::random_vector;
use std::time::Instant;

fn time_dots(v: &Basis<C64>, w: &[C64], k: usize) -> f64 {
    let mut h = vec![C64::ZERO; k];
    let mut best = f64::MAX;
    for _ in 0..20 {
        let t = Instant::now();
        v.view().dots_block(k, w, 1, &mut h, None);
        best = best.min(t.elapsed().as_secs_f64());
    }
    best
}

fn main() {
    let (n, k) = (3000usize, 19usize);
    let mut v = Basis::<C64>::new(n, 2 * k);
    for j in 0..2 * k {
        v.col_mut(j).copy_from_slice(&random_vector::<C64>(n, j as u64 + 1));
    }
    let w = random_vector::<C64>(n, 99);
    let before = time_dots(&v, &w, k);
    let q = random_vector::<C64>(k * (k / 2), 5);
    v.mult_in_place(0, k, k / 2, &q, k); // one GEMM tile call per row tile
    let after = time_dots(&v, &w, k);
    println!(
        "dots_block (n={n}, k={k}, complex): before the GEMM {:.1} us, after {:.1} us (x{:.2})",
        before * 1e6,
        after * 1e6,
        after / before
    );
}
