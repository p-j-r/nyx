// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! S4 study: HODLR compression and the O(n log² n) direct solve on a kernel matrix, and the off-diagonal ranks
//! of the Schur complement on a separator of a 2-D grid (the dense front that a low-rank multifrontal solver
//! would compress). `cargo run --release --example lowrank_study`
use femcore::linalg::Csr;
use femcore::linalg::csr::laplacian_2d;
use femcore::linalg::lowrank::{Hodlr, HodlrOptions};
use femcore::linalg::rng::random_vector;
use femcore::linalg::sparse_lu::{LuOptions, SparseLu};
use std::time::Instant;

fn main() {
    println!("kernel matrix  A_ij = 0.05·ln(|i−j|/n) + diag, tol 1e-8");
    println!(
        "{:>7} {:>10} {:>9} {:>9} {:>9} {:>10} {:>10}",
        "n", "ranks", "storage", "build[s]", "fact[s]", "solve[s]", "residual"
    );
    for n in [1024usize, 2048, 4096, 8192, 16384] {
        let k = move |i: usize, j: usize| {
            if i == j { 2.0 + n as f64 * 0.01 } else { 0.05 * ((i as f64 - j as f64).abs() / n as f64).ln() }
        };
        let t = Instant::now();
        let h = Hodlr::from_entries(n, &k, HodlrOptions { leaf: 64, tol: 1e-8, max_rank: 200 });
        let tb = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let lu = h.factor().unwrap();
        let tf = t.elapsed().as_secs_f64();
        let b = random_vector::<f64>(n, 1);
        let mut x = b.clone();
        let t = Instant::now();
        lu.solve_in_place(&mut x);
        let ts = t.elapsed().as_secs_f64();
        let mut ax = vec![0.0; n];
        h.matvec(&x, &mut ax);
        let res = ax.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt()
            / b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let prof = h.rank_profile();
        println!(
            "{n:>7} {:>10} {:>9.3} {tb:>9.3} {tf:>9.3} {ts:>10.5} {res:>10.1e}",
            prof[0],
            h.storage() as f64 / (n * n) as f64
        );
    }

    println!("\nSchur complement of the 5-point Laplacian on the middle grid line (separator), tol 1e-8");
    println!("{:>6} {:>8} {:>22} {:>9}", "grid", "sep", "ranks per level", "storage");
    for nx in [64usize, 128, 256] {
        let a = laplacian_2d(nx, nx);
        let sc = nx / 2;
        let is_sep = |p: usize| p % nx == sc;
        let interior: Vec<usize> = (0..nx * nx).filter(|&p| !is_sep(p)).collect();
        let sep: Vec<usize> = (0..nx * nx).filter(|&p| is_sep(p)).collect();
        let mut pos = vec![usize::MAX; nx * nx];
        interior.iter().enumerate().for_each(|(k, &p)| pos[p] = k);
        let mut spos = vec![usize::MAX; nx * nx];
        sep.iter().enumerate().for_each(|(k, &p)| spos[p] = k);
        let mut tii: Vec<(usize, usize, f64)> = Vec::new();
        let mut ais: Vec<Vec<(usize, f64)>> = vec![Vec::new(); sep.len()];
        for &p in &interior {
            for (&c, &v) in a.row(p).0.iter().zip(a.row(p).1) {
                if is_sep(c as usize) {
                    ais[spos[c as usize]].push((pos[p], v));
                } else {
                    tii.push((pos[p], pos[c as usize], v));
                }
            }
        }
        let ni = interior.len();
        let lu = SparseLu::factor(&Csr::from_triplets(ni, ni, &tii).unwrap(), &LuOptions::default()).unwrap();
        let ns = sep.len();
        let mut s = vec![0.0; ns * ns];
        for jj in 0..ns {
            let mut col = vec![0.0; ni];
            ais[jj].iter().for_each(|&(i, v)| col[i] = v);
            lu.solve_in_place(&mut col);
            for ii in 0..ns {
                let mut v = if ii == jj {
                    4.0
                } else if ii.abs_diff(jj) == 1 {
                    -1.0
                } else {
                    0.0
                };
                ais[ii].iter().for_each(|&(i, w)| v -= w * col[i]);
                s[jj * ns + ii] = v;
            }
        }
        let h = Hodlr::from_dense(&s, ns, HodlrOptions { leaf: 8, tol: 1e-8, max_rank: ns });
        println!(
            "{nx:>6} {ns:>8} {:>22} {:>9.3}",
            format!("{:?}", h.rank_profile()),
            h.storage() as f64 / (ns * ns) as f64
        );
    }
}
