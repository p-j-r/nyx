// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! SpMV speed and memory: CSR vs block CSR (3×3 blocks, elasticity-like) and CSR vs symmetric
//! (upper-triangle) storage for a complex symmetric matrix. `cargo run --release --example spmv_formats`.
use femcore::linalg::bsr::{Bsr, SymCsr};
use femcore::linalg::csr::laplacian_3d;
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Csr};
use std::time::Instant;

fn blockify(a: &Csr<f64>, b: usize) -> Csr<f64> {
    let mut t = Vec::new();
    for i in 0..a.nrows() {
        let (cs, vs) = a.row(i);
        for (&j, &v) in cs.iter().zip(vs) {
            for r in 0..b {
                for c in 0..b {
                    let w = if r == c { 2.0 + r as f64 } else { 0.25 * (1.0 + (r + c) as f64) };
                    t.push((i * b + r, j as usize * b + c, v * w));
                }
            }
        }
    }
    Csr::from_triplets(a.nrows() * b, a.ncols() * b, &t).unwrap()
}

fn best<F: FnMut()>(mut f: F, reps: usize) -> f64 {
    f();
    (0..reps)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64()
        })
        .fold(f64::MAX, f64::min)
}

fn main() {
    let a3 = laplacian_3d(40);
    for b in [3usize, 2] {
        let a = blockify(&a3, b);
        let bsr = Bsr::from_csr(&a, b).unwrap();
        let x = random_vector::<f64>(a.ncols(), 3);
        let mut y = vec![0.0; a.nrows()];
        let t_csr = best(|| a.spmv(&x, &mut y), 30);
        let t_bsr = best(|| bsr.spmv(&x, &mut y), 30);
        let csr_bytes = a.nnz() * (8 + 4) + (a.nrows() + 1) * 8;
        println!(
            "real, block {b}: n = {}, nnz = {}: CSR {:.2} ms ({:.1} MB) | BSR {:.2} ms ({:.1} MB)  -> {:.2}x faster, {:.2}x smaller",
            a.nrows(),
            a.nnz(),
            t_csr * 1e3,
            csr_bytes as f64 / 1e6,
            t_bsr * 1e3,
            bsr.bytes() as f64 / 1e6,
            t_csr / t_bsr,
            csr_bytes as f64 / bsr.bytes() as f64
        );
    }
    // complex symmetric
    let n = a3.nrows();
    let d: Vec<_> = (0..n).map(|i| (i, i, C64::new(0.0, 0.1 * (i % 5) as f64))).collect();
    let ac = Csr::axpby(C64::ONE, &a3.to_complex(), C64::ONE, &Csr::from_triplets(n, n, &d).unwrap());
    let sym = SymCsr::from_csr(&ac, false, 1e-14).unwrap();
    let x = random_vector::<C64>(n, 3);
    let mut y = vec![C64::ZERO; n];
    let t_csr = best(|| ac.spmv(&x, &mut y), 30);
    let t_sym = best(|| sym.spmv(&x, &mut y), 30);
    let csr_bytes = ac.nnz() * (16 + 4) + (n + 1) * 8;
    println!(
        "complex symmetric 3-D: n = {n}: CSR {:.2} ms ({:.1} MB) | upper triangle {:.2} ms ({:.1} MB) -> {:.2}x time, {:.2}x smaller",
        t_csr * 1e3,
        csr_bytes as f64 / 1e6,
        t_sym * 1e3,
        sym.bytes() as f64 / 1e6,
        t_sym / t_csr,
        csr_bytes as f64 / sym.bytes() as f64
    );
}
