// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! 2-D heat equation `y' + L y = 0` on a 200×200 grid with the variable-step BDF: steps, factorisations
//! and where the time goes (`cargo run --release --example ts_bench`).
use femcore::linalg::csr::laplacian_2d;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::ts::bdf::{BdfOptions, bdf_solve};
use femcore::linalg::ts::newmark::{NewmarkParams, newmark_solve};
use femcore::linalg::{Csr, rng::random_vector};

fn main() {
    if std::env::args().any(|a| a == "reuse") {
        reuse_study();
        return;
    }
    let (nx, ny) = (200, 200);
    let n = nx * ny;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_2d(nx, ny));
    let y0 = random_vector::<f64>(n, 1);
    for rtol in [1e-3, 1e-6] {
        let opts = BdfOptions { definiteness: Definiteness::PositiveDefinite, ..BdfOptions::new(rtol, 1e-8) };
        let sol = bdf_solve(&m, &k, &mut |_, f| f.fill(0.0), 0.0, 10.0, &y0, &[10.0], &opts).unwrap();
        let s = &sol.stats;
        println!(
            "BDF rtol {rtol:e}: {} steps ({} rejected), {} factorisations ({:.3} s), {} solves ({:.3} s), total {:.3} s, steps/order {:?}",
            s.steps,
            s.rejected,
            s.factorizations,
            s.factor_seconds,
            s.solves,
            s.solve_seconds,
            s.seconds,
            s.steps_per_order
        );
    }
    // second-order wave equation, 400 steps, one factorisation
    let sol = newmark_solve(
        &m,
        None,
        &k,
        &mut |_, f| f.fill(0.0),
        0.0,
        40.0,
        400,
        &y0,
        &vec![0.0; n],
        NewmarkParams::generalized_alpha(0.8),
        400,
        Definiteness::PositiveDefinite,
    )
    .unwrap();
    let s = &sol.stats;
    println!(
        "gen-α: {} steps, {} factorisations ({:.3} s), {} solves ({:.3} s), total {:.3} s",
        s.steps, s.factorizations, s.factor_seconds, s.solves, s.solve_seconds, s.seconds
    );
}

/// Reuse of factorisations by refinement vs exact refactorisation, 3-D heat (factorisation ≫ solve).
#[allow(dead_code)]
pub fn reuse_study() {
    let n3 = 40;
    let n = n3 * n3 * n3;
    let (m, k) = (Csr::<f64>::identity(n), femcore::linalg::csr::laplacian_3d(n3));
    let y0 = random_vector::<f64>(n, 1);
    for reuse in [0.0, 0.3] {
        let opts = BdfOptions {
            definiteness: Definiteness::PositiveDefinite,
            reuse_tolerance: reuse,
            ..BdfOptions::new(1e-5, 1e-8)
        };
        let sol = bdf_solve(&m, &k, &mut |_, f| f.fill(0.0), 0.0, 10.0, &y0, &[10.0], &opts).unwrap();
        let s = &sol.stats;
        println!(
            "3-D BDF reuse {reuse}: {} steps, {} factorisations ({:.2} s), total {:.2} s",
            s.steps, s.factorizations, s.factor_seconds, s.seconds
        );
    }
}
