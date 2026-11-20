// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! TOAR vs Krylov–Schur shift-and-invert on the explicit 2n linearisation (time and basis memory).
use femcore::linalg::csr::{laplacian_2d, mass_1d};
use femcore::linalg::direct::Definiteness;
use femcore::linalg::problem::{EigenProblem, KrylovSchurShiftInvert, ShiftInvertEigenSolver, Symmetry};
use femcore::linalg::toar::{QuadraticProblem, ToarOptions, toar_shift_invert};
use femcore::linalg::{C64, Csr};
use std::time::Instant;

fn main() {
    let (nx, ny) = (150usize, 150usize);
    let n = nx * ny;
    let k = laplacian_2d(nx, ny).map(|v| v * 40.0).to_complex();
    let m: Csr<C64> = Csr::identity(n);
    let _ = mass_1d;
    let cm = Csr::axpby(C64::new(0.05, 0.0), &Csr::identity(n), C64::new(0.002, 0.0), &k);
    let p = QuadraticProblem { m, c: cm, k };
    let sigma = C64::new(-0.5, 30.0);
    let t = Instant::now();
    let r = toar_shift_invert(
        &p,
        &ToarOptions::new(6, std::env::var("NCV").ok().and_then(|s| s.parse().ok()).unwrap_or(60), sigma),
    )
    .unwrap();
    println!(
        "TOAR: {:.3} s, converged {}, basis {} scalars (linearised: {}), residuals {:?}",
        t.elapsed().as_secs_f64(),
        r.converged,
        r.basis_scalars,
        r.basis_scalars_linearised,
        r.residuals.iter().fold(0.0f64, |a, &b| a.max(b))
    );
    let (a, b) = p.linearisation();
    let prob = EigenProblem { a, b: Some(b), symmetry: Symmetry::General };
    let t = Instant::now();
    let mut ks = KrylovSchurShiftInvert::new(&prob, Definiteness::General).unwrap();
    ks.opts.tol = 1e-9;
    ks.opts.ncv = Some(60);
    let res = ks.solve_near(sigma, 6).unwrap();
    println!(
        "Krylov–Schur on 2n pencil: {:.3} s, converged {}, matvecs {}, max residual {:?}",
        t.elapsed().as_secs_f64(),
        res.converged,
        res.matvecs,
        res.residuals.iter().fold(0.0f64, |a, &b| a.max(b))
    );
}
